function testMSSDampingConfig()
% Check the sibling MSS loader against saved JSON-backed MAT metadata.
% No installed vessel files are replaced and no BEM solve is run.
root = fileparts(fileparts(mfilename('fullpath')));
originalPath = path;
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() cleanupTest(folder,originalPath)); %#ok<NASGU>
addpath(fullfile(fileparts(root),'MSS','CRAFT','utils_hydroVessel'));
cases = dir(fullfile(root,'vessels_capytaine','*','config.json'));
for i = 1:numel(cases)
    base = cases(i).folder;
    config = jsondecode(fileread(fullfile(base,'config.json')));
    savedPath = fullfile(base,config.output_filename);
    cfg = hydroVesselConfig(savedPath,true);
    expected = config.viscous_damping;
    fields = fieldnames(expected);
    assert(numel(fieldnames(cfg.damping)) == numel(fields));
    for k = 1:numel(fields)
        name = fields{k};
        actualValue = cfg.damping.(name); expectedValue = expected.(name);
        assert(isequal(actualValue(:),expectedValue(:)));
    end
end
% WAMIT vessels now load their importer-selected defaults from the MAT-files.
cfg = hydroVesselConfig('fpso.mat');
assert(isequal(cfg.damping.T_126,[22.7 3.1 3.4]));
assert(isequal(cfg.damping.delta_zeta_345,[0 .1 0]));
assert(isequal(cfg.damping.nonlinear_456,[5 0 0]));
cfg = hydroVesselConfig('semisub.mat');
assert(isequal(cfg.damping.T_126,[6.1 3.1 3.1]));
assert(isequal(cfg.damping.delta_zeta_345,[0 .1 0]));
assert(isequal(cfg.damping.nonlinear_456,[5 0 0]));
% Older custom files without nonlinear metadata default to zero.
vessel.powerBased = struct('T_126',[1 2 3],'delta_zeta_345',[0 .1 0]);
file = fullfile(folder,'custom.mat'); save(file,'vessel','-v7');
cfg = hydroVesselConfig(file,true);
assert(isequal(cfg.damping.nonlinear_456,[0 0 0]));
vessel.powerBased.nonlinear_456 = [5 -1 0]; save(file,'vessel','-v7');
failed = false;
try
    hydroVesselConfig(file,true);
catch problem
    failed = ~isempty(strfind(problem.message,'nonlinear_456'));
end
assert(failed);
fprintf('MSS JSON-backed and WAMIT MAT-file damping defaults passed.\n');
end

function cleanupTest(folder,originalPath)
path(originalPath);
rmdir(folder,'s');
end
