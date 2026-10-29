function testStandalonePeriods()
% Numerical period regression without MSS, graphics, or a BEM solve.
root = fileparts(fileparts(mfilename('fullpath')));
helpers = tempname;
mkdir(helpers);
originalPath = path;
cleanup = onCleanup(@() cleanupHelpers(helpers,originalPath)); %#ok<NASGU>
copyfile(fullfile(root,'matlab','private','computeVesselDamping.m'),helpers);
copyfile(fullfile(root,'matlab','private','computeVesselPeriods.m'),helpers);
addpath(helpers);

base = fullfile(root,'vessels_capytaine','testShip_csv');
config = jsondecode(fileread(fullfile(base,'config.json')));
saved = load(fullfile(base,config.output_filename),'vessel');
vessel = saved.vessel;
vessel.powerBased = config.viscous_damping;
vessel = computeVesselDamping(vessel,0.8);
[T,zeta,omega,omega_n] = computeVesselPeriods(vessel);
% The shared JONSWAP baseline changes the viscous roll damping slightly.
assert(max(abs(T-[3.54720818003 3.52717578396 3.81649045370])) < 1e-6);
assert(max(abs(zeta-[.212395732724 .201724799659 .103266941872])) < 1e-6);
assert(all(omega > 0) && all(omega_n > 0));

% The 10 rad/s endpoint denotes infinity and cannot enter interpolation.
marker = abs(vessel.freqs(:)-10) < 10*eps(10);
assert(nnz(marker)==1);
vessel.A(:,:,marker,1) = 1e20;
vessel.B(:,:,marker,1) = 1e20;
[Tmarker,zetamarker] = computeVesselPeriods(vessel);
assert(max(abs(Tmarker-T)) < 1e-10);
assert(max(abs(zetamarker-zeta)) < 1e-10);
fprintf('Standalone natural periods and infinity-marker exclusion passed.\n');
end

function cleanupHelpers(folder,originalPath)
path(originalPath);
rmdir(folder,'s');
end
