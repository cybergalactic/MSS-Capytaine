function [case_dir, vessel_name] = vesselCatalogueCase(case_name)
% Resolve a grouped catalogue case, accepting group/vessel or a unique name.
case_name = strrep(char(case_name), '\', '/');
if isempty(regexp(case_name, '^[A-Za-z0-9_-]+(/[A-Za-z0-9_-]+)?$', 'once'))
    error('MSSCapytaine:InvalidCase', 'Use a vessel name or group/vessel.');
end
% This helper lives in matlab/private, two folders below the repository root.
root = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'vessels_capytaine');
if any(case_name == '/')
    case_dir = fullfile(root, case_name);
    if exist(fullfile(case_dir, 'config.json'), 'file') ~= 2
        error('MSSCapytaine:UnknownCase', 'Unknown catalogue case: %s', case_name);
    end
else
    matches = dir(fullfile(root, '*', case_name, 'config.json'));
    if numel(matches) ~= 1
        error('MSSCapytaine:UnknownOrAmbiguousCase', ...
            'Use group/vessel to select a unique catalogue case: %s', case_name);
    end
    case_dir = matches(1).folder;
end
[~, vessel_name] = fileparts(case_dir);
end
