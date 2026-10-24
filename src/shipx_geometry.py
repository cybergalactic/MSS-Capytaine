"""Read VERES MGF half sections and clip them to a selected design draft."""

from pathlib import Path
import numpy as np


def read_mgf_sections(path: Path, draft_m: float):
    """Return internal x-aft/y-starboard/z-up sections, waterline origin.

    MGF x points aft from midships; its vertical coordinates are heights
    above baseline. Preserve the contour order, including bulbous sections
    whose vertical coordinates need not be strictly monotone.
    """
    if not np.isfinite(draft_m) or draft_m <= 0:
        raise ValueError("draft_m must be positive and finite for MGF geometry")
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if len(lines) < 8:
        raise ValueError("MGF requires four header lines, Lpp, and station data")
    try:
        lpp = float(lines[4].strip())
        if not np.isfinite(lpp) or lpp <= 0:
            raise ValueError("Invalid MGF Lpp")
        # The first four lines can be blank and must not be stripped away.
        tokens = iter(" ".join(lines[5:]).split())
        sections = []
        while (identifier := next(tokens, None)) is not None:
            int(identifier)
            x = float(next(tokens))
            count = int(next(tokens))
            if count < 2:
                raise ValueError("Each MGF section needs at least two points")
            profile = np.array([[float(next(tokens)), float(next(tokens))]
                                for _ in range(count)])
            if not np.isfinite(x) or not np.isfinite(profile).all() or np.any(profile[:, 0] < 0):
                raise ValueError("MGF sections require finite coordinates and nonnegative half breadth")
            if profile[0, 1] > profile[-1, 1]:
                profile = profile[::-1]
            wet = []
            for p, q in zip(profile[:-1], profile[1:]):
                p_in, q_in = p[1] <= draft_m, q[1] <= draft_m
                if not p_in and not q_in:
                    continue
                if p_in != q_in:
                    crossing = p + (q-p) * (draft_m-p[1]) / (q[1]-p[1])
                    crossing[1] = draft_m
                    p, q = (p, crossing) if p_in else (crossing, q)
                if wet and not np.allclose(wet[-1], p, atol=1e-8, rtol=0):
                    raise ValueError("MGF section has disconnected wetted contours")
                if not wet:
                    wet.append(p)
                if not np.allclose(wet[-1], q, atol=1e-10, rtol=0):
                    wet.append(q)
            if len(wet) < 2:
                continue  # Entirely above water or just tangent to the waterline.
            wet = np.array(wet)
            if np.ptp(wet[:, 1]) <= 1e-10:
                continue
            if not np.isclose(wet[0, 0], 0, atol=1e-6):
                raise ValueError("MGF wetted section must start on the centerline keel")
            if not np.isclose(wet[-1, 1], draft_m, atol=1e-6):
                raise ValueError("MGF section must reach the selected waterline")
            wet[:, 1] -= draft_m
            wet[0, 0] = 0.0
            wet[-1, 1] = 0.0
            sections.append((x, wet))
    except (ValueError, StopIteration) as exc:
        raise ValueError(f"Invalid MGF geometry in {path}: {exc}") from exc
    sections.sort(key=lambda section: section[0])
    if len(sections) < 3 or np.any(np.diff([x for x, _ in sections]) <= 0):
        raise ValueError("MGF needs at least three distinct wetted stations")
    return sections
