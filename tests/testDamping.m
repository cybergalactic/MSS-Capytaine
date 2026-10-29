function testDamping()
% Damping checks use only the MSS-Capytaine helper and saved MAT files.
root = fileparts(fileparts(mfilename('fullpath')));
% MATLAB disallows private folders on the search path. Test an unchanged
% temporary copy instead of exposing the production helper publicly.
helpers = tempname;
mkdir(helpers);
originalPath = path;
cleanup = onCleanup(@() cleanupHelper(helpers,originalPath)); %#ok<NASGU>
copyfile(fullfile(root,'matlab','private','computeVesselDamping.m'),helpers);
addpath(helpers);
cases = {'testShip_csv','otter_55kg_csv', ...
    'otter_80kg_csv','supply_mgf','s175_mgf','tanker_gdf','wigley_stl'};
for c = 1:numel(cases)
    base = fullfile(root,'vessels_capytaine',cases{c});
    config = jsondecode(fileread(fullfile(base,'config.json')));
    saved = load(fullfile(base,config.output_filename));
    v = saved.vessel;
    v.powerBased = config.viscous_damping;
    peaks = [0.8 1.0];
    if ~isempty(strfind(cases{c},'otter_'))
        peaks(end+1) = 1.5; % Regression: current MSS treats explicit inputs as T.
    end
    for wp = peaks
        previousWarning = warning('query','MSS:LegacyKappa126');
        result = computeVesselDamping(v,wp);
        currentWarning = warning('query','MSS:LegacyKappa126');
        assert(strcmp(previousWarning.state,currentWarning.state));
        assert(~isfield(result.powerBased,'kappa_126'));
        m = diag(result.M); d = diag(result.D);
        assert(max(abs(m([1 2 6])./d([1 2 6])-config.viscous_damping.T_126(:))) < 1e-10);
        assert(all(diag(result.powerBased.Bv) >= 0));
        assert(all(result.powerBased.T_126 <= result.powerBased.T_max_126));
        assert(isequal(result.A,v.A) && isequal(result.B,v.B) && isequal(result.C,v.C));
        assert(isequal(result.motionRAO,v.motionRAO));
    end
end
% The submerged call uses the unchanged four-time-constant MSS interface.
base = fullfile(root,'vessels_capytaine','LAUV_marie_csv');
config = jsondecode(fileread(fullfile(base,'config.json')));
auvSaved = load(fullfile(base,config.output_filename));
auv = auvSaved.vessel;
auv.powerBased = config.viscous_damping;
for wp = [0.8 1.5]
    result = computeVesselDamping(auv,wp);
    m = diag(result.M); d = diag(result.D);
    assert(max(abs(m([1 2 3 6])./d([1 2 3 6])-config.viscous_damping.T_1236(:))) < 1e-10);
    assert(all(result.powerBased.T_1236 <= result.powerBased.T_max_1236));
    assert(isequal(result.A,auv.A) && isequal(result.B,auv.B) && isequal(result.C,auv.C));
end
v = saved.vessel;
v.powerBased = struct('T_126',[1 1 1],'delta_zeta_345',[0 .1 0]);
reference = computeVesselDamping(v,.8);
v.powerBased.T_126 = reference.powerBased.T_max_126;
boundary = computeVesselDamping(v,.8);
bv = diag(boundary.powerBased.Bv);
assert(max(abs(bv([1 2 6]))) < 1e-8);
v.B(:) = 0;
v.powerBased.T_126 = [100 200 300];
zero = computeVesselDamping(v,.8);
assert(all(isinf(zero.powerBased.T_max_126)));
v = saved.vessel;
v.powerBased = struct('T_126',[1e6 1e6 1e6],'delta_zeta_345',[0 .1 0]);
try
    computeVesselDamping(v,.8);
    error('Test:ExpectedError','Expected infeasible targets to fail.');
catch problem
    assert(strcmp(problem.identifier,'MSSCapytaine:InfeasibleTimeConstants'));
    assert(~isempty(strfind(problem.message,'surge')));
    assert(~isempty(strfind(problem.message,'sway')));
    assert(~isempty(strfind(problem.message,'yaw')));
    assert(isempty(strfind(problem.message,'omega_p')));
end
v.powerBased = struct('kappa_126',[.05 .05 .05],'delta_zeta_345',[0 .1 0]);
try
    computeVesselDamping(v,.8);
    error('Test:ExpectedError','Expected legacy metadata to fail.');
catch problem
    assert(strcmp(problem.identifier,'MSSCapytaine:LegacyDamping'));
    assert(~isempty(strfind(problem.message,'T_126')));
end
fprintf('Damping targets, maxima, catalogue defaults and migration checks passed.\n');
end

function cleanupHelper(folder,originalPath)
path(originalPath);
rmdir(folder,'s');
end
