results_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    'capytaineTestShip', 'results');
load(fullfile(results_dir, 'capytaineTestShip.mat'));
load(fullfile(results_dir, 'generated_hull_panels.mat'));

omega_p = 0.8;
vessel = computeManeuveringModel(vessel,omega_p,0);

% vesselPeriods expects frequency-indexed A/B and one 6-by-6 C matrix.
Aw = vessel.A(:,:,:,1);
Bw = vessel.B(:,:,:,1) + vessel.Bv(:,:,:,1);
C = vessel.C(:,:,1,1);

[T,zeta,omega,omega_n] = vesselPeriods( ...
    vessel.freqs, vessel.MRB, Aw, Bw, C, 'coupled', true);

figure(gcf)
patch('Vertices', vertices, ...
    'Faces', faces, ...
    'FaceColor', [0.8 0.8 0.8], ...
    'EdgeColor', 'k');

axis equal
grid on
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')
view(3)

plotTF(vessel,'force','rads')
plotABC(vessel,'A')
plotABC(vessel,'B')
plotBv(vessel)
