function geometry = readVesselGeometry(case_name)
% Read source CSV/MGF/GDF/STL without a solve. Returned coordinates are MSS FSD.
% Accept a catalogue case (group/vessel or unique name), or a JSON file path.
% isfile checks an actual path; exist(...,'file') also finds MSS functions
% such as tanker.m or supply.m on the MATLAB path for a short vessel name.
if isfile(char(case_name))
    config_path = char(case_name);
    case_dir = fileparts(config_path);
    if isempty(case_dir), case_dir = pwd; end
    [~, vessel_name] = fileparts(case_dir);
else
    [case_dir, vessel_name] = vesselCatalogueCase(case_name);
    config_path = fullfile(case_dir, 'config.json');
end
config = jsondecode(fileread(config_path));
if isfield(config, 'body_name'), vessel_name = char(config.body_name); end
if isfield(config, 'geometry_file')
    filename = char(config.geometry_file);
elseif isfield(config, 'offset_points_csv')
    filename = char(config.offset_points_csv);
else
    error('MSSCapytaine:MissingGeometry', 'Expected geometry_file in config.json.');
end
path = fullfile(case_dir, filename);
if exist(filename, 'file') == 2 && isAbsolute(filename), path = filename; end
if exist(path, 'file') ~= 2
    error('MSSCapytaine:MissingGeometry', 'Geometry file not found: %s', path);
end
format = detectFormat(path);
geometry = struct('name', vessel_name, 'format', format, 'path', path, ...
    'config_path', config_path, 'results_dir', fullfile(case_dir,'results'), ...
    'vertices', [], 'faces', [], 'sections', {{}}, 'wetted_sections', {{}}, ...
    'center_of_mass_m', [], 'section_centers_m', 0);
if isfield(config,'demihull_spacing_m')
    spacing=config.demihull_spacing_m;
    if ~strcmp(format,'csv')
        error('MSSCapytaine:InvalidSpacing','demihull_spacing_m is only supported for CSV demi-hull offsets.');
    end
    if ~isnumeric(spacing) || ~isreal(spacing) || ~isscalar(spacing) || ...
            ~isfinite(spacing) || spacing<=0
        error('MSSCapytaine:InvalidSpacing','demihull_spacing_m must be a positive finite centerline distance.');
    end
    geometry.section_centers_m=[-spacing/2,spacing/2];
end
if isfield(config,'output_dir')
    geometry.results_dir=fullfile(case_dir,char(config.output_dir));
    if isAbsolute(char(config.output_dir)), geometry.results_dir=char(config.output_dir); end
end
if isfield(config,'center_of_mass_m')
    cg=config.center_of_mass_m;
    if numel(cg)~=3 || any(~isfinite(cg(:)))
        error('MSSCapytaine:InvalidCG','center_of_mass_m must contain three finite FSD coordinates.');
    end
    geometry.center_of_mass_m=cg(:)';
end
switch format
    case 'csv'
        fid = fopen(path, 'rt'); cleanup = onCleanup(@() fclose(fid));
        header = strsplit(stripBom(fgetl(fid)), ',');
        header = cellfun(@strtrim, header, 'UniformOutput', false);
        names = {'x_m','z_m','half_breadth_m'}; indices = zeros(1,3);
        for i = 1:3
            index = find(strcmp(header, names{i}));
            if numel(index) ~= 1, error('MSSCapytaine:InvalidCSV', 'Missing column %s', names{i}); end
            indices(i) = index;
        end
        values = textscan(fid, repmat('%f',1,numel(header)), 'Delimiter', ',', 'CollectOutput', true);
        values = values{1};
        if isempty(values) || any(~isfinite(values(:)))
            error('MSSCapytaine:InvalidCSV', 'CSV coordinates must be finite.');
        end
        x = values(:,indices(1)); z = values(:,indices(2)); y = values(:,indices(3));
        if any(y < 0), error('MSSCapytaine:InvalidCSV', 'Half breadth must be nonnegative.'); end
        stations = unique(x);
        if numel(stations) < 3, error('MSSCapytaine:InvalidCSV', 'At least three stations are required.'); end
        submerged = isfield(config,'submerged') && config.submerged;
        for i = 1:numel(stations)
            rows = x == stations(i); p = [x(rows),y(rows),z(rows)];
            if size(p,1) < 2 || any(diff(p(:,3)) > 1e-9) || abs(p(1,2)) > 1e-9
                error('MSSCapytaine:InvalidCSV', 'Sections must run from the centerline bottom to top (decreasing FSD z).');
            end
            if (~submerged && (any(p(:,3)<-1e-9) || abs(p(end,3))>1e-9)) || ...
                    (submerged && abs(p(end,2))>1e-9)
                error('MSSCapytaine:InvalidCSV', 'Invalid section endpoint or FSD vertical coordinates.');
            end
            geometry.sections{end+1} = p;
        end
        geometry.wetted_sections = geometry.sections;
        if isfield(config,'demihull_spacing_m') && config.demihull_spacing_m<=2*max(y)
            error('MSSCapytaine:InvalidSpacing','demihull_spacing_m must exceed the demi-hull beam; hulls must not touch or overlap.');
        end
    case 'mgf'
        if ~isfield(config,'draft_m') || ~isscalar(config.draft_m) || ...
                ~isfinite(config.draft_m) || config.draft_m <= 0
            error('MSSCapytaine:InvalidDraft', 'MGF geometry requires positive draft_m.');
        end
        if isfield(config,'submerged') && config.submerged
            error('MSSCapytaine:InvalidMGF', 'MGF import currently supports surface vessels only.');
        end
        draft = config.draft_m;
        fid = fopen(path,'rt'); cleanup = onCleanup(@() fclose(fid));
        for i=1:4, if ~ischar(fgetl(fid)), error('MSSCapytaine:InvalidMGF','Incomplete MGF header.'); end, end
        lpp = fscanf(fid,'%f',1);
        if isempty(lpp) || ~isfinite(lpp) || lpp<=0, error('MSSCapytaine:InvalidMGF','Invalid MGF Lpp.'); end
        while true
            identifier = fscanf(fid,'%f',1);
            if isempty(identifier), break; end
            x = fscanf(fid,'%f',1); count = fscanf(fid,'%f',1);
            if isempty(x) || isempty(count) || ~isfinite(x) || ~isfinite(count) || count<2 || count~=fix(count)
                error('MSSCapytaine:InvalidMGF','Invalid station header.');
            end
            profile = fscanf(fid,'%f',[2,count])';
            if size(profile,1)~=count || any(~isfinite(profile(:))) || any(profile(:,1)<0)
                error('MSSCapytaine:InvalidMGF','Invalid section coordinates.');
            end
            if profile(1,2)>profile(end,2), profile=flipud(profile); end
            geometry.sections{end+1} = [repmat(-x,count,1),profile(:,1),draft-profile(:,2)];
            wet = clipProfile(profile, draft);
            if ~isempty(wet)
                geometry.wetted_sections{end+1} = [repmat(-x,size(wet,1),1),wet(:,1),draft-wet(:,2)];
            end
        end
        if numel(geometry.wetted_sections)<3
            error('MSSCapytaine:InvalidMGF','MGF needs at least three wetted stations.');
        end
    case 'stl'
        fid=fopen(path,'rb','ieee-le'); cleanup=onCleanup(@() fclose(fid));
        fseek(fid,0,'eof'); bytes=ftell(fid); fseek(fid,80,'bof'); count=fread(fid,1,'uint32');
        if ~isempty(count) && count>0 && bytes==84+50*double(count)
            vertices=zeros(3*count,3);
            for i=1:count
                row=fread(fid,12,'single=>double'); fread(fid,1,'uint16');
                vertices(3*i-2:3*i,:)=reshape(row(4:12),3,3)';
            end
        else
            text=fileread(path);
            facets=regexp(text,'facet\s+normal([\s\S]*?)endfacet','tokens','ignorecase');
            if isempty(facets), error('MSSCapytaine:InvalidSTL','No STL facets found.'); end
            vertices=zeros(3*numel(facets),3);
            for i=1:numel(facets)
                points=regexp(facets{i}{1},'vertex\s+([^\r\n]+)','tokens','ignorecase');
                if numel(points)~=3, error('MSSCapytaine:InvalidSTL','Each facet needs three vertices.'); end
                for j=1:3
                    p=sscanf(points{j}{1},'%f');
                    if numel(p)~=3, error('MSSCapytaine:InvalidSTL','Each vertex needs three coordinates.'); end
                    vertices(3*i-3+j,:)=p';
                end
            end
        end
        if any(~isfinite(vertices(:))), error('MSSCapytaine:InvalidSTL','STL coordinates must be finite.'); end
        geometry.vertices=vertices; geometry.faces=reshape(1:size(vertices,1),3,[])';
    case 'gdf'
        fid=fopen(path,'rt'); cleanup=onCleanup(@() fclose(fid)); fgetl(fid);
        scales=sscanf(fgetl(fid),'%f'); flags=sscanf(fgetl(fid),'%f'); count=sscanf(fgetl(fid),'%f',1);
        if numel(scales)<2 || any(scales(1:2)<=0) || numel(flags)<2 || ...
                any(~ismember(flags(1:2),[0 1])) || isempty(count) || count<1 || count~=fix(count)
            error('MSSCapytaine:InvalidGDF','Invalid GDF header.');
        end
        values=fscanf(fid,'%f');
        if numel(values)~=12*count || any(~isfinite(values))
            error('MSSCapytaine:InvalidGDF','GDF panel count does not match its coordinates.');
        end
        geometry.vertices=reshape(values,3,[])' * diag([1 -1 -1]);
        geometry.faces=reshape(1:4*count,4,[])';
        for dimension=1:2
            if flags(dimension)==1
                reflected=geometry.vertices; reflected(:,dimension)=-reflected(:,dimension);
                faces=geometry.faces(:,end:-1:1)+size(geometry.vertices,1);
                geometry.vertices=[geometry.vertices;reflected]; geometry.faces=[geometry.faces;faces];
            end
        end
end
end

function result=isAbsolute(path)
result=~isempty(path) && (path(1)=='/' || ~isempty(regexp(path,'^[A-Za-z]:[\\/]','once')));
end

function format=detectFormat(path)
[~,~,extension]=fileparts(path); extension=lower(extension);
if any(strcmp(extension,{'.csv','.mgf','.gdf','.stl'})), format=extension(2:end); return; end
fid=fopen(path,'rb','ieee-le'); cleanup=onCleanup(@() fclose(fid));
fseek(fid,0,'eof'); bytes=ftell(fid); fseek(fid,80,'bof'); count=fread(fid,1,'uint32');
if ~isempty(count) && count>0 && bytes==84+50*double(count), format='stl'; return; end
frewind(fid); prefix=char(fread(fid,4096,'*uint8')');
lines=regexp(stripBom(prefix),'\r\n|\n|\r','split');
if all(ismember({'x_m','z_m','half_breadth_m'},strtrim(strsplit(lines{1},','))))
    format='csv'; return
end
if ~isempty(regexpi(prefix,'facet\s+normal')) && ~isempty(regexpi(prefix,'vertex\s+')) && ~isempty(regexpi(prefix,'endfacet'))
    format='stl'; return
end
if numel(lines)>=4
    scales=sscanf(lines{2},'%f'); flags=sscanf(lines{3},'%f'); count=sscanf(lines{4},'%f');
    if numel(scales)>=2 && all(scales(1:2)>0) && numel(flags)>=2 && all(ismember(flags(1:2),[0 1])) && numel(count)==1 && count>0 && count==fix(count)
        format='gdf'; return
    end
end
if numel(lines)>=9
    lpp=sscanf(lines{5},'%f'); identifier=sscanf(lines{6},'%f'); x=sscanf(lines{7},'%f'); count=sscanf(lines{8},'%f'); point=sscanf(lines{9},'%f');
    if numel(lpp)==1 && lpp>0 && numel(identifier)==1 && identifier>=0 && identifier==fix(identifier) && numel(x)==1 && numel(count)==1 && count>=2 && count==fix(count) && numel(point)==2
        format='mgf'; return
    end
end
error('MSSCapytaine:UnknownGeometry','Cannot identify geometry contents; use CSV, MGF, GDF or STL.');
end

function wet=clipProfile(profile,draft)
wet=[];
for i=1:size(profile,1)-1
    p=profile(i,:); q=profile(i+1,:); pin=p(2)<=draft; qin=q(2)<=draft;
    if ~pin && ~qin, continue; end
    if pin~=qin
        crossing=p+(q-p)*(draft-p(2))/(q(2)-p(2)); crossing(2)=draft;
        if pin, q=crossing; else, p=crossing; end
    end
    if ~isempty(wet) && norm(wet(end,:)-p,inf)>1e-8
        error('MSSCapytaine:InvalidMGF','Disconnected wetted section contour.');
    end
    if isempty(wet), wet=p; end
    if norm(wet(end,:)-q,inf)>1e-10, wet=[wet;q]; end
end
if size(wet,1)<2 || max(wet(:,2))-min(wet(:,2))<=1e-10, wet=[]; return; end
if abs(wet(1,1))>1e-6 || abs(wet(end,2)-draft)>1e-6
    error('MSSCapytaine:InvalidMGF','Wetted contour must run from centerline keel to waterline.');
end
wet(1,1)=0; wet(end,2)=draft;
end

function text=stripBom(text)
if numel(text)>=3 && isequal(double(text(1:3)),[239 187 191])
    text=text(4:end);
elseif ~isempty(text) && double(text(1))==65279
    text=text(2:end);
end
end
