function testDemihullGeometry()
% Check local CSV sections and configured twin-hull placement without graphics.
original_path=path;
restore_path=onCleanup(@() path(original_path));
addpath(fullfile(fileparts(which('plotVesselGeometry')),'private'));
a=readVesselGeometry('otter_80kg_csv');
b=readVesselGeometry('otter_80kg');
assert(strcmp(a.path,b.path)); assert(strcmp(a.format,'csv'));
assert(isequal(a.section_centers_m,[-0.395 0.395]));
points=vertcat(a.sections{:});
assert(abs(max(points(:,2))-0.125)<1e-9);
assert(abs(max(points(:,3))-0.1951219512)<1e-9);
assert(abs(a.center_of_mass_m(3)+0.10386139215)<5.1e-5);
light=readVesselGeometry('otter_55kg');
assert(isequal(light.section_centers_m,a.section_centers_m));
assert(abs(light.center_of_mass_m(1)-0.2010045508)<5.1e-5);
assert(abs(light.center_of_mass_m(3)+0.1989634709)<5.1e-5);
assert(max(points(:,2))+a.section_centers_m(1)<0);
assert(-max(points(:,2))+a.section_centers_m(2)>0);
% A normal CSV remains a single local centerline.
single=readVesselGeometry('testShip'); assert(isequal(single.section_centers_m,0));
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s'));
config=jsondecode(fileread(a.config_path)); config.geometry_file=a.path;
file=fullfile(folder,'config.json');
for spacing={0,-1,true,[],0.25,0.1}
    config.demihull_spacing_m=spacing{1};
    fid=fopen(file,'wt'); fprintf(fid,'%s',jsonencode(config)); fclose(fid);
    failed=false;
    try, readVesselGeometry(file);
    catch err, failed=strcmp(err.identifier,'MSSCapytaine:InvalidSpacing'); end
    assert(failed);
end
stl=readVesselGeometry('wigley');
config.geometry_file=stl.path; config.demihull_spacing_m=0.79;
fid=fopen(file,'wt'); fprintf(fid,'%s',jsonencode(config)); fclose(fid);
failed=false;
try, readVesselGeometry(file);
catch err, failed=strcmp(err.identifier,'MSSCapytaine:InvalidSpacing'); end
assert(failed);
disp('CSV demi-hulls: short-name lookup, local sections, spacing and invalid-input checks passed.');
end
