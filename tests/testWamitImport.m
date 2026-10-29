function testWamitImport()
% Import all bundled WAMIT cases in a temporary folder; no solver is run.
suiteRoot = fileparts(fileparts(mfilename('fullpath')));
root = fullfile(fileparts(suiteRoot),'MSS','HYDRO','vessels_wamit');
addpath(root);
oldDirectory = pwd;
oldPath = path;
pauseState = pause('query');
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() restoreTest(oldDirectory,oldPath,pauseState,folder)); %#ok<NASGU>
pause('off');
importerFolder = fullfile(folder,'importer'); mkdir(importerFolder);
copyfile(fullfile(root,'wamit2vessel.m'),fullfile(importerFolder,'wamit2vessel.m'));
addpath(importerFolder,'-begin');
assert(~isempty(strfind(which('wamit2vessel'), ...
    [filesep 'importer' filesep 'wamit2vessel.m'])));
names = {'tanker','fpso','semisub'};
dimensions = [10 246 46; 12 200 44; 21 115 80];
timeConstants = [49.6 3.6 4.5; 22.7 3.1 3.4; 6.1 3.1 3.1];
extensions = {'.1','.3','.4','.9','.out','.frc','.pot'};
for i = 1:numel(names)
    name = names{i}; source = fullfile(root,[name '_wamitData']);
    work = fullfile(folder,name); mkdir(work);
    for j = 1:numel(extensions)
        file = [name extensions{j}];
        copyfile(fullfile(source,file),fullfile(work,file));
    end
    cd(folder);
    vessel = wamit2vessel(fullfile(name,name), ...
        dimensions(i,1),dimensions(i,2),dimensions(i,3), ...
        '0000',timeConstants(i,:),[0 .1 0],[5 0 0]);
    outputFile = fullfile(importerFolder,[name '.mat']);
    assert(exist(outputFile,'file')==2);
    assert(exist(fullfile(work,[name '.mat']),'file')~=2);
    saved = load(outputFile);
    assert(isequal(saved.vessel,vessel));
    assert(strcmp(vessel.main.name,name));
    assert(~isfield(vessel,'driftfrc') && isfield(vessel,'drift'));
    assert(numel(vessel.velocities)==1 && vessel.velocities==0);
    assert(all(isfinite(vessel.A(:))) && all(isfinite(vessel.B(:))));
    zeroFrequencyDamping = vessel.B(:,:,1);
    assert(all(zeroFrequencyDamping(:)==0));
    reference = load(fullfile(root,[name '.mat']));
    assert(isequal(vessel.drift,reference.vessel.drift));
    assert(isequal(vessel.powerBased,reference.vessel.powerBased));
    assert(~isfield(vessel.powerBased,'A_eq') && ~isfield(vessel.powerBased,'Bv'));
    for k = 1:6
        assert(size(vessel.forceRAO.amp{k},3)==1);
        assert(size(vessel.motionRAO.amp{k},3)==1);
    end
    fprintf('Full WAMIT import passed: %s.\n',name);
end
end

function restoreTest(oldDirectory,oldPath,pauseState,folder)
cd(oldDirectory);
path(oldPath);
pause(pauseState);
rmdir(folder,'s');
end
