function vessel = computeVesselDamping(vessel, omega_p)
% Compute power-based equivalent matrices and apply catalogue time constants.
% This uses only the saved vessel structure and MATLAB/Octave functions.

if nargin < 2 || isempty(omega_p)
    omega_p = 1.0;
end

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
G = zeros(6);
G([3 4 5],[3 4 5]) = vessel.C([3 4 5],[3 4 5],1);
heaveIsRestored = abs(G(3,3)) > 1e-10 * max(1,norm(G,'fro'));
if submerged == heaveIsRestored
    error('MSSCapytaine:InvalidRestoring', ...
        'The saved submerged setting and heave restoring matrix disagree.');
end
vessel = powerBasedEquivalent(vessel,omega_p,G);
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
for k = 1:numel(restored)
    i = restored(k);
    vessel.powerBased.Bv(i,i) = 2 * ratios(k) * sqrt(m(i) * G(i,i));
end
vessel.D = vessel.powerBased.B_eq + vessel.powerBased.Bv;
vessel.powerBased.(timeField) = T;
vessel.powerBased.(ratioField) = ratios;
vessel.powerBased.(['T_max_' timeField(3:end)]) = Tmax;
end

function vessel = powerBasedEquivalent(vessel,omega_p,G)
% The spectrum scale cancels when it is normalized; only its shape matters.
if ~isnumeric(omega_p) || ~isreal(omega_p) || ~isscalar(omega_p) || ...
        ~isfinite(omega_p) || omega_p <= 0
    error('MSSCapytaine:InvalidPeakFrequency', ...
        'omega_p must be a positive finite scalar in rad/s.');
end
w = vessel.freqs(:);
finite = isfinite(w) & w < 10;
w = w(finite);
if numel(w) < 2 || any(w < 0) || any(diff(w) <= 0)
    error('MSSCapytaine:InvalidFrequencies', ...
        'At least two increasing finite frequencies below the infinity marker are required.');
end
% Match the JONSWAP baseline in MSS computeManeuveringModel.
wFine = linspace(max(w(1),1e-6),w(end),1000)';
gamma = 3.3;
sigma = 0.09 * ones(size(wFine));
sigma(wFine < omega_p) = 0.07;
peak = exp(-0.5 * ((wFine - omega_p) ./ (sigma * omega_p)).^2);
weight = wFine.^(-5) .* exp(-1.25 * (omega_p ./ wFine).^4) .* gamma.^peak;
m0 = trapz(wFine,weight);
if ~isfinite(m0) || m0 <= 0
    error('MSSCapytaine:InvalidSpectrum', ...
        'The wave-spectrum weight is zero or nonfinite on this frequency grid.');
end
weight = weight / m0;
Aeq = zeros(6); Beq = zeros(6);
for i = 1:6
    for j = 1:6
        Aij = reshape(vessel.A(i,j,finite,1),[],1);
        Bij = reshape(vessel.B(i,j,finite,1),[],1);
        Aeq(i,j) = trapz(wFine,interp1(w,Aij,wFine,'pchip') .* weight);
        Beq(i,j) = trapz(wFine,interp1(w,Bij,wFine,'pchip') .* weight);
    end
end
vessel.powerBased.omega_p = omega_p;
vessel.powerBased.spectrumType = 'JONSWAP';
vessel.powerBased.gamma = gamma;
vessel.powerBased.A_eq = Aeq;
vessel.powerBased.B_eq = Beq;
vessel.powerBased.Bv = zeros(6);
vessel.MA = Aeq;
vessel.M = vessel.MRB + Aeq;
if rank(vessel.M) < 6
    error('MSSCapytaine:SingularInertia','The equivalent inertia matrix is singular.');
end
[Q,R,P] = qr(vessel.M);
vessel.Minv = P * (R \ Q');
vessel.G = G;
vessel.D = Beq;
end
