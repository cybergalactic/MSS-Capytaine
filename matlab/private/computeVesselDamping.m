function vessel = computeVesselDamping(vessel, omega_p)
% Apply catalogue time constants using MSS equivalent matrices.
% The saved-input call below supports both the older percentage interface
% and the current time-constant interface, without ambiguous explicit inputs.

p = vessel.powerBased;
if isfield(p, 'kappa_126')
    error('MSSCapytaine:LegacyDamping', ...
        ['Legacy kappa_126 percentage settings are no longer supported. ' ...
         'Replace them in config.json with T_126 = [T_surge T_sway T_yaw] ' ...
         'in seconds, and run python main.py <case> --update-damping. ' ...
         'Do not copy [0.05 0.05 0.05] as time constants. No BEM rerun is needed.']);
end
submerged = isfield(vessel.main, 'submerged') && logical(vessel.main.submerged);
if submerged
    timeField = 'T_1236'; ratioField = 'delta_zeta_45';
    idx = [1 2 3 6]; restored = [4 5];
else
    timeField = 'T_126'; ratioField = 'delta_zeta_345';
    idx = [1 2 6]; restored = [3 4 5];
end
if ~isfield(p,timeField) || ~isfield(p,ratioField)
    error('MSSCapytaine:InvalidDamping', ...
        'Saved damping inputs require %s and %s. Refresh them from config.json.', timeField,ratioField);
end
T = p.(timeField); ratios = p.(ratioField);
if ~isnumeric(T) || ~isreal(T) || ~isvector(T) || numel(T) ~= numel(idx) || ...
        any(~isfinite(T(:))) || any(T(:) <= 0)
    error('MSSCapytaine:InvalidDamping','%s must contain %d positive finite time constants in seconds.',timeField,numel(idx));
end
if ~isnumeric(ratios) || ~isreal(ratios) || ~isvector(ratios) || ...
        numel(ratios) ~= numel(restored) || any(~isfinite(ratios(:))) || any(ratios(:) < 0)
    error('MSSCapytaine:InvalidDamping','%s must contain %d nonnegative finite damping-ratio increments.',ratioField,numel(restored));
end
T = reshape(T,1,[]); ratios = reshape(ratios,1,[]);
if submerged
    % Keep the established submerged formulation in the installed MSS.
    vessel = computeManeuveringModel(vessel,omega_p,T,ratios,0);
else
    % Ask for radiation alone through the saved-percentage compatibility
    % interface, which both MSS versions support. An explicit [0 0 0] is
    % invalid with the current MSS because explicit inputs now mean seconds.
    % Apply and validate the user's actual T_126 below, reporting all failures
    % together. The internal zero increments are not user damping settings.
    vessel.powerBased.kappa_126 = [0 0 0];
    vessel.powerBased.delta_zeta_345 = ratios;
    previousWarning = warning('query','MSS:LegacyKappa126');
    warning('off','MSS:LegacyKappa126');
    restoreWarning = onCleanup(@() warning(previousWarning));
    vessel = computeManeuveringModel(vessel,omega_p);
    clear restoreWarning
    vessel.powerBased = rmfield(vessel.powerBased,'kappa_126');
end
m = diag(vessel.M); b = diag(vessel.powerBased.B_eq);
if any(~isfinite(m(idx))) || any(m(idx) <= 0) || ...
        any(~isfinite(b(idx))) || any(b(idx) < 0)
    error('MSSCapytaine:InvalidEquivalentModel', ...
        'Equivalent inertia must be positive and radiation damping nonnegative in the unrestrained motions.');
end
Tmax = inf(size(T)); positive = b(idx)' > 0;
Tmax(positive) = m(idx(positive))' ./ b(idx(positive))';
target = m(idx)' ./ T;
extra = target - b(idx)';
if any(~isfinite(target))
    error('MSSCapytaine:InvalidDamping','Time constants are too small to produce finite damping coefficients.');
end
tolerance = 1e-10 * max([ones(numel(idx),1),abs(target(:)),abs(b(idx))],[],2)';
bad = find(extra < -tolerance);
if ~isempty(bad)
    names = {'surge','sway','heave','roll','pitch','yaw'};
    message = '';
    for index = 1:numel(bad)
        k = bad(index);
        message = [message sprintf(['Requested %s time constant: %.6g s\n' ...
            'Maximum supported by the potential-flow equivalent model: %.6g s\n\n'], ...
            names{idx(k)},T(k),Tmax(k))]; %#ok<AGROW>
    end
    error('MSSCapytaine:InfeasibleTimeConstants','%s%s',message, ...
        ['The requested time constants would require negative additional viscous damping ' ...
         'because the computed radiation damping alone produces a faster response. ' ...
         'Specify time constants no greater than the corresponding maxima, or review ' ...
         'the model and operating conditions if the values come from measurements.']);
end
for k = 1:numel(idx)
    vessel.powerBased.Bv(idx(k),idx(k)) = max(0,extra(k));
end
vessel.D = vessel.powerBased.B_eq + vessel.powerBased.Bv;
vessel.powerBased.(timeField) = T;
vessel.powerBased.(ratioField) = ratios;
vessel.powerBased.(['T_max_' timeField(3:end)]) = Tmax;
end
