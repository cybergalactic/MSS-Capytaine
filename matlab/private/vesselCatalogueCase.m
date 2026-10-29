function [case_dir, vessel_name] = vesselCatalogueCase(case_name)
% Resolve a vessel-format folder or a unique JSON body_name.
case_name = strrep(char(case_name), '\', '/');
if isempty(regexp(case_name, '^[A-Za-z0-9_-]+$', 'once'))
    error('MSSCapytaine:InvalidCase', 'Use a vessel name or vessel-format folder.');
end
% This helper lives in matlab/private, two folders below the repository root.
root = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'vessels_capytaine');
case_dir = fullfile(root, case_name);
if exist(fullfile(case_dir, 'config.json'), 'file') == 2
    config = jsondecode(fileread(fullfile(case_dir, 'config.json')));
else
    entries = dir(fullfile(root, '*', 'config.json'));
    matches = {};
    for i = 1:numel(entries)
        candidate = jsondecode(fileread(fullfile(entries(i).folder, 'config.json')));
        if isfield(candidate, 'body_name') && strcmp(char(candidate.body_name), case_name)
            matches{end+1} = entries(i).folder;
        end
    end
    if isempty(matches)
        error('MSSCapytaine:UnknownCase', 'Unknown catalogue case: %s', case_name);
    elseif numel(matches) ~= 1
        error('MSSCapytaine:UnknownOrAmbiguousCase', ...
            'Use a vessel-format folder to select a unique catalogue case: %s', case_name);
    end
    case_dir = matches{1};
    config = jsondecode(fileread(fullfile(case_dir, 'config.json')));
end
vessel_name = char(config.body_name);
end
