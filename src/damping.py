"""Validate JSON damping settings and refresh saved MAT metadata without BEM.

T_126/T_1236 are total diagonal velocity-decay targets in seconds. Python
checks their syntax; MATLAB checks feasibility against the equivalent radiation
damping, which depends on the selected post-processing spectrum.
"""

from pathlib import Path
import json
import os
import tempfile

import numpy as np
from scipy.io import loadmat, savemat


def validate_damping(settings: object, submerged: bool = False) -> dict[str, np.ndarray]:
    """Return validated time constants and restored-mode ratio increments."""
    if not isinstance(settings, dict):
        raise ValueError("viscous_damping must be an object")
    if "kappa_126" in settings:
        raise ValueError(
            "Legacy viscous_damping.kappa_126 is no longer supported. "
            "Replace the percentage array (for example [0.05, 0.05, 0.05]) "
            "with T_126: [T_surge, T_sway, T_yaw] in seconds. "
            "Choose positive total velocity-decay time constants; do not "
            "copy the percentage values. Keep delta_zeta_345 for heave, "
            "roll and pitch. Then run python main.py <case> --update-damping. "
            "No BEM rerun is needed for damping-only changes."
        )
    schema = ({"T_1236": (4, True), "delta_zeta_45": (2, False)} if submerged
              else {"T_126": (3, True), "delta_zeta_345": (3, False)})
    unknown, missing = set(settings) - set(schema), set(schema) - set(settings)
    if unknown:
        raise ValueError("Unknown viscous_damping entries: " + ", ".join(sorted(unknown)))
    if missing:
        raise ValueError("Missing viscous_damping entries: " + ", ".join(sorted(missing)))
    parameters = {}
    for name, (size, positive) in schema.items():
        raw = np.asarray(settings[name], dtype=object)
        if raw.shape != (size,) or any(isinstance(x, (bool, str)) for x in raw.flat):
            raise ValueError(f"viscous_damping.{name} must contain {size} finite numbers")
        try:
            values = np.asarray(settings[name], dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError(f"viscous_damping.{name} must contain {size} finite numbers") from error
        if not np.all(np.isfinite(values)):
            raise ValueError(f"viscous_damping.{name} must contain {size} finite numbers")
        if positive and np.any(values <= 0):
            raise ValueError(f"viscous_damping.{name} must be positive")
        if not positive and np.any(values < 0):
            raise ValueError(f"viscous_damping.{name} must be nonnegative")
        parameters[name] = values
    return parameters


def refresh_damping(config_path: Path) -> Path:
    """Replace only a saved vessel's damping metadata from its JSON config.

    Preserve raw MATLAB structures and cells, including all hydrodynamic arrays.
    This does not rerun BEM or certify changes to geometry/mass/frequency inputs.
    """
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(config.get("submerged"), bool):
        raise ValueError("submerged must be true or false")
    parameters = validate_damping(config.get("viscous_damping"), config["submerged"])
    result = config_path.parent / config.get("output_dir", "results") / config["output_filename"]
    if not result.is_file():
        raise FileNotFoundError(f"No saved vessel file to update: {result}. Run the BEM case first.")
    data = {k: v for k, v in loadmat(result).items() if not k.startswith("__")}
    vessel = data.get("vessel")
    if vessel is None or vessel.shape != (1, 1) or "powerBased" not in (vessel.dtype.names or ()):
        raise ValueError("The saved MAT file must contain a scalar vessel with powerBased metadata")
    power = np.empty((1, 1), dtype=[(key, object) for key in parameters])
    for key, value in parameters.items():
        power[key][0, 0] = value.reshape(1, -1)
    vessel["powerBased"][0, 0] = power
    # A processed file can contain derived caches. Never leave stale results
    # after changing inputs; the raw A/B/C/MRB and RAOs remain untouched.
    cached = {"M", "MA", "Minv", "D", "G"}
    names = [key for key in vessel.dtype.names if key not in cached]
    if len(names) != len(vessel.dtype.names):
        fresh = np.empty((1, 1), dtype=[(key, object) for key in names])
        for key in names:
            fresh[key][0, 0] = vessel[key][0, 0]
        data["vessel"] = fresh
    fd, temporary = tempfile.mkstemp(dir=result.parent, suffix=".mat")
    os.close(fd)
    try:
        savemat(temporary, data, long_field_names=True, do_compression=True)
        os.chmod(temporary, result.stat().st_mode & 0o777)
        os.replace(temporary, result)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return result
