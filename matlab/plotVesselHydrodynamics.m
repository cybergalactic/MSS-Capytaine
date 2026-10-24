function vessel = plotVesselHydrodynamics(vessel_name, omega_p)
% plotVesselHydrodynamics processes and plots MSS-Capytaine hydrodynamic data.
%
% The vessel structure is loaded automatically from:
%
%   vessels_capytaine/<format>/myVessel/results/myVessel.mat
%
% The optional omega_p argument selects the MSS power-based weighting peak
% in rad/s (default 0.8, an illustrative choice rather than a vessel constant).
% Choose the physical sea-state peak, not a radiation-damping resonance.
% It uses the damping inputs exported from config.json.
% T_126/T_1236 specify total diagonal velocity-decay time constants in seconds.
% Infeasible targets produce an error with the potential-flow maximum values.
% After editing only damping settings, use python main.py <case> --update-damping.
% It displays the constant mass, damping, and restoring matrices and plots the
% force RAOs, added mass, radiation damping, and viscous damping. Natural
% periods and damping ratios are also reported for a surface vessel.
%
% Examples:
%   vessel = plotVesselHydrodynamics('testShip');
%
%   vessel = plotVesselHydrodynamics('LAUV_marie');
%   vessel = plotVesselHydrodynamics('otter_55kg', 1.5);
% Changing omega_p reweights the saved A/B data; no Python rerun is needed.
% This post-processing choice does not add forward-speed hydrodynamics.
%
% MSS dependencies:
%   computeManeuveringModel.m
%   vesselPeriods.m
%   plotTF.m
%   plotABC.m
%   plotBv.m
%
% Author: Thor I. Fossen
% Date: 2026-10-07

[case_dir, vessel_name] = vesselCatalogueCase(vessel_name);
if nargin < 2 || isempty(omega_p), omega_p = 0.8; end
if ~isnumeric(omega_p) || ~isreal(omega_p) || ~isscalar(omega_p) || ...
        ~isfinite(omega_p) || omega_p <= 0
    error('MSSCapytaine:InvalidPeakFrequency','omega_p must be a positive finite scalar in rad/s.');
end
vessel_file = fullfile(case_dir, 'results', [vessel_name '.mat']);

if exist(vessel_file, 'file') ~= 2
    error('MSSCapytaine:VesselFileNotFound', ...
        ['Could not find the generated vessel file:\n%s\n' ...
         'Run "python main.py %s" from the MSS-Capytaine repository root.'], ...
        vessel_file, vessel_name);
end

vessel_data = load(vessel_file, 'vessel');
if ~isfield(vessel_data, 'vessel') || ~isstruct(vessel_data.vessel)
    error('MSSCapytaine:InvalidVesselFile', ...
        'The generated MATLAB file must contain an MSS vessel structure.');
end
vessel = vessel_data.vessel;

required_fields = {'main', 'MRB', 'A', 'B', 'C', 'freqs', ...
    'forceRAO', 'powerBased'};
for index = 1:numel(required_fields)
    if ~isfield(vessel, required_fields{index})
        error('MSSCapytaine:InvalidVessel', ...
            'The vessel structure is missing field "%s".', ...
            required_fields{index});
    end
end

if ~isfield(vessel.main, 'name') || ...
        ~strcmp(char(vessel.main.name), vessel_name)
    error('MSSCapytaine:InvalidVesselFile', ...
        'The loaded vessel.main.name must be "%s".', vessel_name);
end

vessel = computeVesselDamping(vessel, omega_p);

if isfield(vessel.powerBased,'T_1236')
    dampingDofs = [1 2 3 6];
    targetT = vessel.powerBased.T_1236;
    maximumT = vessel.powerBased.T_max_1236;
else
    dampingDofs = [1 2 6];
    targetT = vessel.powerBased.T_126;
    maximumT = vessel.powerBased.T_max_126;
end
names = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
fprintf('\nMotion       Specified T [s]     T_max [s]\n');
for index = 1:numel(dampingDofs)
    fprintf('%-10s %14.6g %14.6g\n',names{dampingDofs(index)},targetT(index),maximumT(index));
end

display(vessel.M, 'M')
display(vessel.D, 'D')
display(vessel.G, 'G')

is_submerged = isfield(vessel.main, 'submerged') && ...
    logical(vessel.main.submerged);
if ~is_submerged
    vesselPeriods(vessel.freqs, vessel.MRB, vessel.A, ...
        vessel.B + vessel.powerBased.Bv, vessel.G, 'coupled', true);
end

plotTF(vessel, 'force', 'rads')
plotABC(vessel, 'A')
plotABC(vessel, 'B')
plotBv(vessel)

end
