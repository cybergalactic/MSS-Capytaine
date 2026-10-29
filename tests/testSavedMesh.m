function testSavedMesh()
% Reader checks without graphics or a BEM solve.
original_path=path; restore_path=onCleanup(@() path(original_path));
addpath(fullfile(fileparts(which('plotVesselGeometry')),'private'));
folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s'));
fid=fopen(fullfile(folder,'hull.stl'),'wt');
fprintf(fid,['solid hull\nfacet normal 0 0 1\nouter loop\n' ...
    'vertex 0 0 0\nvertex 1 0 0\nvertex 0 1 0\nendloop\nendfacet\nendsolid hull\n']); fclose(fid);
config=struct('geometry_file','hull.stl','body_name','example', ...
    'output_filename','different_name.mat', ...
    'output_dir','custom_results','center_of_mass_m',[0 0 9]);
file=fullfile(folder,'config.json'); writeConfig(file,config);
geometry=readVesselGeometry(file);
assert(strcmp(geometry.results_dir,fullfile(folder,'custom_results')));
mesh=readCapytaineMesh(geometry); assert(~mesh.available);
mkdir(geometry.results_dir);
vertices=[0 0 0;1 0 0;1 1 1;0 1 1]; faces=[1 2 3 4;1 2 3 3];
coordinate_system='MSS_FSD';
save(mesh.path,'vertices','faces','coordinate_system');
vessel=struct('main',struct('CG',[0 0 0.25]));
save(fullfile(folder,'different_name.mat'),'vessel');
mesh=readCapytaineMesh(geometry);
assert(mesh.available && isequal(mesh.vertices,vertices) && isequal(mesh.faces,faces));
assert(isequal(mesh.center_of_mass_m,[0 0 0.25])); % Saved CG, not edited JSON CG.
% Absolute output folders work too.
config.output_dir=geometry.results_dir; writeConfig(file,config);
absolute=readVesselGeometry(file); assert(strcmp(absolute.results_dir,geometry.results_dir));
% Legacy solver-axis files are rotated to FSD, preserving face connectivity.
save(mesh.path,'vertices','faces'); legacy=readCapytaineMesh(geometry);
assert(isequal(legacy.vertices,vertices*diag([-1 1 -1])));
assert(isequal(legacy.faces,faces));
coordinate_system='unknown'; save(mesh.path,'vertices','faces','coordinate_system');
assertError(@() readCapytaineMesh(geometry),'MSSCapytaine:InvalidAxes');
coordinate_system='MSS_FSD'; faces(1,1)=5;
save(mesh.path,'vertices','faces','coordinate_system');
assertError(@() readCapytaineMesh(geometry),'MSSCapytaine:InvalidFaces');
save(mesh.path,'vertices');
assertError(@() readCapytaineMesh(geometry),'MSSCapytaine:InvalidMeshFile');
disp('Saved mesh: missing results, custom folders, saved CG, FSD/legacy axes and invalid panels passed.');
end

function writeConfig(file,config)
fid=fopen(file,'wt'); fprintf(fid,'%s',jsonencode(config)); fclose(fid);
end

function assertError(action,identifier)
failed=false;
try, action(); catch err, failed=strcmp(err.identifier,identifier); end
assert(failed);
end
