function mesh = readCapytaineMesh(geometry)
% Read saved solver hull panels in FSD; a missing solve is a normal preview.
mesh = struct('available',false, ...
    'path',fullfile(geometry.results_dir,'generated_hull_panels.mat'), ...
    'vertices',[],'faces',[],'center_of_mass_m',[]);
if ~isfile(mesh.path), return; end
data=load(mesh.path,'vertices','faces','coordinate_system');
if ~isfield(data,'vertices') || ~isfield(data,'faces')
    error('MSSCapytaine:InvalidMeshFile','The saved mesh must contain vertices and faces.');
end
vertices=data.vertices; faces=data.faces;
if ~isnumeric(vertices) || size(vertices,2)~=3 || isempty(vertices) || ...
        any(~isfinite(vertices(:)))
    error('MSSCapytaine:InvalidVertices','Hull-panel vertices must be a finite numeric N-by-3 array.');
end
if ~isnumeric(faces) || ~ismember(size(faces,2),[3 4]) || isempty(faces) || ...
        any(~isfinite(faces(:))) || any(faces(:)~=fix(faces(:))) || ...
        any(faces(:)<1) || any(faces(:)>size(vertices,1))
    error('MSSCapytaine:InvalidFaces','Hull panels must contain valid MATLAB triangle or quadrilateral vertex indices.');
end
if ~isfield(data,'coordinate_system')
    % Older exports stored internal x-aft, y-starboard, z-up coordinates.
    vertices=vertices*diag([-1 1 -1]);
elseif ~strcmp(char(data.coordinate_system),'MSS_FSD')
    error('MSSCapytaine:InvalidAxes','Expected MSS_FSD panel coordinates.');
end
mesh.vertices=vertices; mesh.faces=double(faces); mesh.available=true;
% Use the solved CG rather than a potentially edited JSON CG on the mesh side.
vessel_path=geometry.vessel_path;
if isfile(vessel_path)
    data=load(vessel_path,'vessel');
    if isfield(data,'vessel') && isfield(data.vessel,'main') && isfield(data.vessel.main,'CG')
        cg=data.vessel.main.CG;
        if isnumeric(cg) && numel(cg)==3 && all(isfinite(cg(:)))
            mesh.center_of_mass_m=cg(:)';
        end
    end
end
end
