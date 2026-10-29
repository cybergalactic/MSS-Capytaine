function [T,zeta,omega,omega_n] = computeVesselPeriods(vessel)
% Coupled heave, roll and pitch periods from saved frequency-dependent data.
% The artificial 10 rad/s infinity marker is excluded from interpolation.

w = vessel.freqs(:);
finite = isfinite(w) & abs(w-10) >= 10*eps(10);
w = w(finite);
if numel(w) < 2 || any(w < 0) || any(diff(w) <= 0)
    error('MSSCapytaine:InvalidFrequencies', ...
        'Periods require two increasing nonnegative finite frequencies.');
end
nf = numel(w);
A = reshape(vessel.A(:,:,finite,1),36,nf).';
B = reshape(vessel.B(:,:,finite,1),36,nf).';
MRB = vessel.MRB;
G = vessel.G;
Bv = vessel.powerBased.Bv;
dofs = [3 4 5];
T = zeros(1,3); zeta = T; omega = T; omega_n = T;

for k = 1:3
    target = dofs(k);
    constraint = @(om) om^2 - selectedEigenvalue(om,target,w,A,MRB,G);
    omega_n(k) = fzero(constraint,[w(1) w(end)]);
    M = MRB + interpolateMatrix(w,A,omega_n(k));
    radiation = interpolateMatrix(w,B,omega_n(k));
    [V,D] = eig(G,M);
    lambda = real(diag(D));
    mode = selectMode(V,lambda,target);
    x = V(:,mode) / max(norm(V(:,mode)),eps);
    modalMass = real(x' * M * x);
    modalStiffness = real(x' * G * x);
    modalDamping = real(x' * (radiation+Bv) * x);
    zeta(k) = modalDamping / (2*sqrt(modalStiffness*modalMass));
    omega(k) = omega_n(k) * sqrt(max(0,1-zeta(k)^2));
    T(k) = 2*pi / omega(k);
end
end

function lambda = selectedEigenvalue(om,target,w,A,MRB,G)
M = MRB + interpolateMatrix(w,A,om);
[V,D] = eig(G,M);
values = real(diag(D));
lambda = values(selectMode(V,values,target));
end

function mode = selectMode(V,lambda,target)
valid = find(isfinite(lambda) & lambda > 1e-9);
if isempty(valid)
    error('MSSCapytaine:NoRestoredMode', ...
        'No positive finite eigenvalue was found for a restored mode.');
end
[~,index] = max(abs(V(target,valid)));
mode = valid(index);
end

function matrix = interpolateMatrix(w,columns,omega)
matrix = reshape(interp1(w,columns,omega,'pchip','extrap'),6,6);
end
