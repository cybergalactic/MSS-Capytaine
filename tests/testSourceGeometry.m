function geometries = testSourceGeometry()
% Run from MATLAB/Octave with addpath('matlab'); addpath('tests'). No graphics.
% Check private-helper access through the public functions before exposing
% internals temporarily for the reader tests below.
entry_points={@plotVesselGeometry,@plotVesselHydrodynamics};
for i=1:numel(entry_points)
    failed=false;
    try, entry_points{i}('unknown_format/unknown_vessel');
    catch err, failed=strcmp(err.identifier,'MSSCapytaine:UnknownCase'); end
    assert(failed);
end
original_path=path;
restore_path=onCleanup(@() path(original_path));
addpath(fullfile(fileparts(which('plotVesselGeometry')),'private'));
cases={'capytaine_offsets_csv/testShip','capytaine_offsets_csv/LAUV_marie', ...
    'shipx_mgf/supply','wamit_gdf/tanker','mesh_stl/wigley'};
expected={'csv','csv','mgf','gdf','stl'};
geometries=cell(size(cases));
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s'));
% MSS provides tanker.m and supply.m; these must not be treated as JSON files.
for name={'tanker','supply'}
    fid=fopen(fullfile(folder,[name{1} '.m']),'wt');
    fprintf(fid,'function value = %s()\nvalue = [];\nend\n',name{1}); fclose(fid);
end
addpath(folder);
assert(exist('tanker','file')==2 && ~isfile('tanker'));
assert(exist('supply','file')==2 && ~isfile('supply'));
for i=1:numel(cases)
    a=readVesselGeometry(cases{i}); assert(strcmp(a.format,expected{i}));
    [~,short_name]=fileparts(cases{i});
    short=readVesselGeometry(short_name);
    assert(strcmp(short.path,a.path));
    assert(isequal(short.vertices,a.vertices)); assert(isequal(short.sections,a.sections));
    assert(~isempty(a.vertices) || ~isempty(a.sections));
    mesh=readCapytaineMesh(a); assert(mesh.available);
    assert(size(mesh.vertices,2)==3 && ~isempty(mesh.faces));
    [case_dir,~]=vesselCatalogueCase(cases{i});
    config=jsondecode(fileread(fullfile(case_dir,'config.json')));
    target=fullfile(folder,'geometry_without_extension'); copyfile(a.path,target);
    config.geometry_file='geometry_without_extension';
    file=fullfile(folder,'config.json'); fid=fopen(file,'wt'); fprintf(fid,'%s',jsonencode(config)); fclose(fid);
    b=readVesselGeometry(file); assert(strcmp(b.format,a.format));
    assert(isequal(a.vertices,b.vertices)); assert(isequal(a.faces,b.faces));
    assert(isequal(a.sections,b.sections)); assert(isequal(a.wetted_sections,b.wetted_sections));
    missing_mesh=readCapytaineMesh(b); assert(~missing_mesh.available);
    geometries{i}=a;
end
% Also test ASCII STL without an extension (the catalogue case is binary).
fid=fopen(target,'wt'); fprintf(fid,['solid hull\nfacet normal 0 0 1\nouter loop\n' ...
    'vertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid hull\n']); fclose(fid);
a=readVesselGeometry(file); assert(strcmp(a.format,'stl')); assert(size(a.faces,1)==1);
fid=fopen(target,'wt'); fprintf(fid,'Unrecognized geometry'); fclose(fid);
failed=false;
try, readVesselGeometry(file); catch err, failed=strcmp(err.identifier,'MSSCapytaine:UnknownGeometry'); end
assert(failed);
disp('Source geometry: four formats, short names with MSS function collisions, missing extensions, ASCII/binary STL passed.');
end
