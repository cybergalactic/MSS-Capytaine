function testWamitMeanDrift()
% Check current saved drift coefficients against the original WAMIT .9 rows.
suiteRoot = fileparts(fileparts(mfilename('fullpath')));
root = fullfile(fileparts(suiteRoot),'MSS','HYDRO','vessels_wamit');
addpath(root);
names = {'tanker','fpso','semisub'};
for i = 1:numel(names)
    folder = fullfile(root,[names{i} '_wamitData']);
    source = fullfile(folder,[names{i} '.9']);
    data = load(fullfile(root,[names{i} '.mat']));
    vessel = data.vessel;
    drift = readWamitMeanDrift(source,vessel.main.rho,vessel.main.g,1);
    assert(isfield(vessel,'drift') && ~isfield(vessel,'driftfrc'));
    for k = 1:6
        assert(size(vessel.forceRAO.amp{k},3)==numel(vessel.velocities));
        assert(size(vessel.motionRAO.amp{k},3)==numel(vessel.velocities));
    end
    assert(isequal(drift.w,vessel.drift.w));
    raw = load(source,'-ascii');
    dofs = [1 2 6]; signs = [1 1 -1]; parity = [1 -1 -1];
    for j = 1:3
        values = vessel.drift.coeff{j};
        assert(size(values,1)==numel(drift.w) && size(values,2)==36);
        assert(size(values,3)==numel(vessel.velocities));
        assert(isequal(values(:,:,1),drift.coeff{j}));
        assert(all(isfinite(values(:))));
        assert(isequal(values(:,20:36),parity(j)*values(:,18:-1:2)));
        rows = raw(raw(:,4)==dofs(j),:);
        for k = 1:size(rows,1)
            [distance,freq] = min(abs(drift.w-2*pi/rows(k,1)));
            assert(distance<1e-12);
            heading = 1+rows(k,2)/10;
            expected = signs(j)*double(vessel.main.rho)*double(vessel.main.g)*rows(k,7);
            assert(abs(values(freq,heading)-expected)<=1e-12*max(1,abs(expected)));
        end
    end
    assert(strcmp(vessel.drift.reference,'WAMIT body origin'));
end
% Synthetic data exercise non-unit length scaling, signs and rejected input.
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
file = fullfile(folder,'synthetic.8');
rows = [];
for period = [2 4]
    for heading = 0:10:180
        for dof = [1 2 6]
            rows(end+1,:) = [period heading heading dof 2 180 -2 0]; %#ok<AGROW>
        end
    end
end
dlmwrite(file,rows,'precision',17);
drift = readWamitMeanDrift(file,1000,10,3);
assert(drift.coeff{1}(1,1)==-60000);
assert(drift.coeff{2}(1,1)==-60000);
assert(drift.coeff{3}(1,1)==180000);
bad = rows; bad(1,3) = 10;
dlmwrite(file,bad,'precision',17); assertRejected(file);
dlmwrite(file,rows(2:end,:),'precision',17); assertRejected(file);
dlmwrite(file,[rows; rows(1,:)],'precision',17); assertRejected(file);
% Old drift fields must never be plotted by the current plotter.
legacy.driftfrc = struct('w',1,'amp',{{1,2,3}});
failed = false;
try
    plotWD(legacy);
catch problem
    failed = ~isempty(strfind(problem.message,'vessel.drift'));
end
assert(failed);
fprintf('WAMIT mean-drift import, saved data, scaling and validation passed.\n');
end

function assertRejected(file)
failed = false;
try
    readWamitMeanDrift(file,1000,10,3);
catch problem
    failed = strcmp(problem.identifier,'MSS:WamitMeanDrift');
end
assert(failed);
end

