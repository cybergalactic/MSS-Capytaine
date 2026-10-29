function vessel = plotVesselHydrodynamics(vessel_name, omega_p)
% plotVesselHydrodynamics processes and plots MSS-Capytaine hydrodynamic data.
%
% The vessel structure is loaded automatically from:
%
%   vessels_capytaine/<vessel_format>/<output_filename>
%
% The optional omega_p argument selects the MSS power-based weighting peak
% in rad/s (default 1.0, the fixed MSS reference rather than a vessel constant).
% An explicit value can be used to inspect another reference sea state.
% It uses the damping inputs exported from config.json.
% T_126/T_1236 specify total diagonal velocity-decay time constants in seconds.
% Infeasible targets produce an error with the potential-flow maximum values.
% After editing only damping settings, use python main.py <case> --update-damping.
% It displays the constant mass, damping, and restoring matrices and plots the
% force RAOs, signed mean wave-drift coefficients when present, added mass,
% radiation damping, and viscous damping, in that order. Natural
% periods and damping ratios are also reported for a surface vessel.
% MATLAB groups the plots into two tabbed windows: wave loads and coefficients.
%
% Examples:
%   vessel = plotVesselHydrodynamics('testShip');
%
%   vessel = plotVesselHydrodynamics('LAUV_marie');
%   vessel = plotVesselHydrodynamics('otter_55kg', 1.5);
% Changing omega_p reweights the saved A/B data; no Python rerun is needed.
% This post-processing choice does not add forward-speed hydrodynamics.
%
% The viewer uses only this repository's saved MAT file and local helpers;
% no MSS installation is needed.
%
% Author: Thor I. Fossen
% Date: 2026-10-07

[case_dir, vessel_name] = vesselCatalogueCase(vessel_name);
if nargin < 2 || isempty(omega_p), omega_p = 1.0; end
if ~isnumeric(omega_p) || ~isreal(omega_p) || ~isscalar(omega_p) || ...
        ~isfinite(omega_p) || omega_p <= 0
    error('MSSCapytaine:InvalidPeakFrequency','omega_p must be a positive finite scalar in rad/s.');
end
config = jsondecode(fileread(fullfile(case_dir, 'config.json')));
vessel_file = fullfile(case_dir, char(config.output_filename));

if exist(vessel_file, 'file') ~= 2
    error('MSSCapytaine:VesselFileNotFound', ...
        ['Could not find the generated vessel file:\n%s\n' ...
         'Run "python main.py %s" from the MSS-Capytaine repository root.'], ...
        vessel_file, vessel_name);
end

vessel_data = load(vessel_file, 'vessel');
if ~isfield(vessel_data, 'vessel') || ~isstruct(vessel_data.vessel)
    error('MSSCapytaine:InvalidVesselFile', ...
        'The generated MATLAB file must contain an MSS vessel structure.');
end
vessel = vessel_data.vessel;

required_fields = {'main', 'MRB', 'A', 'B', 'C', 'freqs', ...
    'forceRAO', 'powerBased'};
for index = 1:numel(required_fields)
    if ~isfield(vessel, required_fields{index})
        error('MSSCapytaine:InvalidVessel', ...
            'The vessel structure is missing field "%s".', ...
            required_fields{index});
    end
end

if ~isfield(vessel.main, 'name') || ...
        ~strcmp(char(vessel.main.name), vessel_name)
    error('MSSCapytaine:InvalidVesselFile', ...
        'The loaded vessel.main.name must be "%s".', vessel_name);
end

vessel = computeVesselDamping(vessel, omega_p);

if isfield(vessel.powerBased,'T_1236')
    dampingDofs = [1 2 3 6];
    targetT = vessel.powerBased.T_1236;
    maximumT = vessel.powerBased.T_max_1236;
else
    dampingDofs = [1 2 6];
    targetT = vessel.powerBased.T_126;
    maximumT = vessel.powerBased.T_max_126;
end
names = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
fprintf('\nMotion       Specified T [s]     T_max [s]\n');
for index = 1:numel(dampingDofs)
    fprintf('%-10s %14.6g %14.6g\n',names{dampingDofs(index)},targetT(index),maximumT(index));
end

display(vessel.M, 'M')
display(vessel.D, 'D')
display(vessel.G, 'G')

is_submerged = isfield(vessel.main, 'submerged') && ...
    logical(vessel.main.submerged);
if ~is_submerged
    try
        [periods,zeta,omega,omega_n] = computeVesselPeriods(vessel);
        fprintf('\n    DOF      T [s]      zeta     omega [rad/s]  omega_n [rad/s]\n');
        fprintf('    ---   ----------  --------  -------------   ---------------\n');
        for index = 1:3
            fprintf('  %4d  %10.3f %10.3f   %10.3f     %10.3f\n', ...
                index+2,periods(index),zeta(index),omega(index),omega_n(index));
        end
    catch problem
        reason = lower(problem.message);
        if isempty(strfind(reason,'bracket')) && isempty(strfind(reason,'sign'))
            rethrow(problem);
        end
        warning('MSSCapytaine:PeriodsUnavailable', ...
            'Could not determine coupled natural periods: %s',problem.message);
    end
end

plotHydrodynamicData(vessel,vessel_name);

end

function plotHydrodynamicData(vessel,vessel_name)
% Draw directly into two tabbed MATLAB figures. Octave uses named figures.
tabbed = ~exist('OCTAVE_VERSION','builtin') && usejava('desktop');
if tabbed
    [waveFigure,waveTabs] = tabWindow('MSSCapytaineWaveLoads', ...
        sprintf('%s: wave loads',vessel_name),[80 80 1180 780]);
    [coefficientFigure,coefficientTabs] = tabWindow('MSSCapytaineCoefficients', ...
        sprintf('%s: coefficients',vessel_name),[130 40 1280 860]);
else
    waveTabs = []; coefficientTabs = [];
end

names = {'Surge','Sway','Heave','Roll','Pitch','Yaw'};
for dof = 1:6
    label = [names{dof} ' RAO'];
    parent = plotPage(waveTabs,label,['MSSCapytaineForce' num2str(dof)], ...
        [vessel_name ': ' label],tabbed);
    plotForceRAO(parent,vessel,dof,names{dof});
end
if isfield(vessel,'drift')
    parent = plotPage(waveTabs,'Wave drift','MSSCapytaineDrift', ...
        [vessel_name ': wave drift'],tabbed);
    plotMeanDrift(parent,vessel);
end

labels = {'Added mass','Potential damping'};
letters = {'A','B'};
for kind = 1:2
    for speed = 1:numel(vessel.velocities)
        label = labels{kind};
        if numel(vessel.velocities) > 1
            label = sprintf('%s %.3g',label,vessel.velocities(speed));
        end
        parent = plotPage(coefficientTabs,label, ...
            sprintf('MSSCapytaine%s%d',letters{kind},speed), ...
            [vessel_name ': ' label],tabbed);
        plotCoefficientMatrix(parent,vessel,letters{kind},speed,tabbed);
    end
end
parent = plotPage(coefficientTabs,'Viscous lateral damping', ...
    'MSSCapytaineDampingLateral',[vessel_name ': lateral damping'],tabbed);
plotEquivalentDamping(parent,vessel,[2 4 6]);
parent = plotPage(coefficientTabs,'Viscous longitudinal damping', ...
    'MSSCapytaineDampingLongitudinal',[vessel_name ': longitudinal damping'],tabbed);
plotEquivalentDamping(parent,vessel,[1 3 5]);

if tabbed
    set(coefficientFigure,'Visible','on');
    set(waveFigure,'Visible','on');
    figure(waveFigure);
end
end

function [fig,tabs] = tabWindow(tag,name,position)
fig = findall(0,'Type','figure','Tag',tag);
if isempty(fig)
    fig = figure('IntegerHandle','off','Tag',tag,'Visible','off', ...
        'WindowStyle','normal','NumberTitle','off','Name',name, ...
        'Position',position);
else
    fig = fig(1);
    set(fig,'Visible','off','WindowStyle','normal','Name',name, ...
        'NumberTitle','off');
    clf(fig);
end
tabs = uitabgroup(fig,'Units','normalized','Position',[0 0 1 1]);
end

function parent = plotPage(tabs,label,tag,name,tabbed)
if tabbed
    parent = uitab(tabs,'Title',label);
else
    parent = findall(0,'Type','figure','Tag',tag);
    if isempty(parent)
        parent = figure('Name',name,'NumberTitle','off','Tag',tag);
    else
        parent = parent(1);
        clf(parent);
        set(parent,'Name',name,'NumberTitle','off');
        figure(parent);
    end
end
end

function plotForceRAO(parent,vessel,dof,name)
w = vessel.forceRAO.w(:);
headings = mod(vessel.headings(:)*180/pi,360);
wanted = 0:30:180;
chosen = zeros(size(wanted));
for k = 1:numel(wanted)
    index = find(abs(headings-wanted(k)) < 1e-6,1);
    if isempty(index)
        error('MSSCapytaine:MissingHeading', ...
            'Force RAO plot requires the %g-degree heading.',wanted(k));
    end
    chosen(k) = index;
end
amplitude = vessel.forceRAO.amp{dof}(:,:,1);
phase = vessel.forceRAO.phase{dof}(:,:,1);
if size(amplitude,1) ~= numel(w) || size(phase,1) ~= numel(w)
    error('MSSCapytaine:InvalidForceRAO', ...
        'Force RAO frequency and coefficient dimensions disagree.');
end
labels = arrayfun(@(beta) sprintf('%g deg',beta),wanted,'UniformOutput',false);
unit = 'kN';
if dof > 3, unit = 'kNm'; end
top = axes('Parent',parent,'Units','normalized','Position',[.09 .56 .86 .36]);
plot(top,w,amplitude(:,chosen)/1000,'LineWidth',1.5);
title(top,['Force RAO amplitude: ' upper(name)]);
ylabel(top,unit); grid(top,'on');
legend(top,labels,'Location','best','FontSize',9);
bottom = axes('Parent',parent,'Units','normalized','Position',[.09 .10 .86 .34]);
plot(bottom,w,(180/pi)*unwrap(phase(:,chosen),[],1),'LineWidth',1.5);
title(bottom,'Phase'); xlabel(bottom,'Wave frequency (rad/s)');
ylabel(bottom,'deg'); grid(bottom,'on');
end

function plotMeanDrift(parent,vessel)
drift = vessel.drift;
if ~isstruct(drift) || ~isfield(drift,'w') || ~isfield(drift,'coeff') || ...
        ~iscell(drift.coeff) || numel(drift.coeff) ~= 3
    error('MSSCapytaine:InvalidDrift', ...
        'Expected vessel.drift.w and three vessel.drift.coeff cells.');
end
w = drift.w(:);
headings = vessel.headings(:);
if ~isnumeric(w) || ~isreal(w) || isempty(w) || any(~isfinite(w)) || ...
        any(w <= 0) || any(w >= 10) || any(diff(w) <= 0) || ...
        ~isnumeric(headings) || ~isreal(headings) || isempty(headings) || ...
        any(~isfinite(headings))
    error('MSSCapytaine:InvalidDrift', ...
        'Mean drift requires positive finite frequencies and finite headings.');
end
for j = 1:3
    values = drift.coeff{j};
    if ~isnumeric(values) || ~isreal(values) || any(~isfinite(values(:))) || ...
            size(values,1) ~= numel(w) || size(values,2) ~= numel(headings) || ...
            size(values,3) < 1 || ndims(values) > 3
        error('MSSCapytaine:InvalidDrift', ...
            'Drift cells must contain real finite frequency-by-heading-by-speed arrays.');
    end
end
degrees = headings*180/pi;
colors = jet(numel(headings));
labels = arrayfun(@(beta) sprintf('%g deg',beta),degrees,'UniformOutput',false);
names = {'SURGE','SWAY','YAW'};
units = {'N/m^2','N/m^2','N m/m^2'};
for j = 1:3
    ax = axes('Parent',parent,'Units','normalized', ...
        'Position',[.1 .72-.25*(j-1) .86 .18]);
    lines = plot(ax,w,drift.coeff{j}(:,:,1),'LineWidth',1.1);
    for k = 1:numel(headings)
        style = '-';
        if degrees(k) >= 180-1e-8, style = '--'; end
        set(lines(k),'Color',colors(k,:),'LineStyle',style);
    end
    title(ax,['Signed mean wave-drift coefficient: ' names{j}]);
    xlabel(ax,'Wave frequency (rad/s)'); ylabel(ax,units{j});
    if numel(w) == 1
        xlim(ax,w(1)*[.95 1.05]);
    else
        xlim(ax,[w(1) w(end)]);
    end
    grid(ax,'on'); set(ax,'FontSize',11);
end
leg = legend(ax,lines,labels,'Location','southoutside', ...
    'NumColumns',min(9,numel(headings)),'FontSize',9);
set(leg,'Units','normalized','Position',[.1 .02 .86 .12]);
end

function plotCoefficientMatrix(parent,vessel,letter,speed,tabbed)
w = vessel.freqs(:);
infinite = abs(w-10) < 10*eps(10);
physical = isfinite(w) & ~infinite;
H = vessel.(letter);
groups = {[1 3 5],[2 4 6]};
for block = 1:2
    indices = groups{block};
    for row = 1:3
        for column = 1:3
            i = indices(row); j = indices(column);
            displayRow = row + 4*(block-1);
            ax = matrixAxes(parent,displayRow,column);
            values = reshape(H(i,j,:,speed),[],1);
            plot(ax,w(physical),values(physical),'b-o','MarkerSize',3);
            if any(infinite)
                hold(ax,'on');
                plot(ax,w(infinite),values(infinite),'rx', ...
                    'MarkerSize',7,'LineWidth',1.5);
                hold(ax,'off');
            end
            title(ax,sprintf('%s_{%d%d}',letter,i,j),'Interpreter','tex');
            if row == 3, xlabel(ax,'rad/s'); end
            grid(ax,'on'); set(ax,'FontSize',8);
        end
    end
end
if tabbed
    annotation(parent,'textbox',[.15 .94 .7 .035], ...
        'String','Longitudinal: surge, heave, pitch','EdgeColor','none', ...
        'HorizontalAlignment','center','FontWeight','bold');
    annotation(parent,'textbox',[.15 .48 .7 .035], ...
        'String','Lateral: sway, roll, yaw','EdgeColor','none', ...
        'HorizontalAlignment','center','FontWeight','bold');
end
end

function ax = matrixAxes(parent,row,column)
x = [.07 .38 .69];
y = [.81 .69 .57 NaN .36 .24 .12];
ax = axes('Parent',parent,'Units','normalized', ...
    'Position',[x(column) y(row) .27 .085]);
end

function plotEquivalentDamping(parent,vessel,indices)
w = vessel.freqs(:);
infinite = abs(w-10) < 10*eps(10);
physical = isfinite(w) & ~infinite;
Beq = vessel.powerBased.B_eq;
Bv = vessel.powerBased.Bv;
x = [.07 .38 .69]; y = [.68 .42 .16];
for row = 1:3
    for column = 1:3
        i = indices(row); j = indices(column);
        ax = axes('Parent',parent,'Units','normalized', ...
            'Position',[x(column) y(row) .25 .20]);
        values = reshape(vessel.B(i,j,:,1),[],1);
        handles = plot(ax,w(physical),values(physical),'b-o','MarkerSize',3);
        hold(ax,'on');
        labels = {'B potential'};
        if any(infinite)
            handles(end+1) = plot(ax,w(infinite),values(infinite),'ro', ...
                'MarkerFaceColor','b','MarkerSize',6);
            labels{end+1} = 'B(\infty)';
        end
        handles(end+1) = plot(ax,[0 10],Beq(i,j)*[1 1],'c-.','LineWidth',1.5);
        labels{end+1} = 'B_{eq}';
        handles(end+1) = plot(ax,[0 10],Bv(i,j)*[1 1],'r-','LineWidth',1.5);
        labels{end+1} = 'B_v';
        handles(end+1) = plot(ax,[0 10],(Beq(i,j)+Bv(i,j))*[1 1], ...
            'm-.','LineWidth',1.5);
        labels{end+1} = 'D total';
        hold(ax,'off');
        title(ax,sprintf('B_{%d%d}',i,j),'Interpreter','tex');
        if row == 3, xlabel(ax,'rad/s'); end
        grid(ax,'on'); set(ax,'FontSize',9);
        if row == 1 && column == 1
            legendAxes = ax;
            legendHandles = handles;
            legendLabels = labels;
        end
    end
end
leg = legend(legendAxes,legendHandles,legendLabels, ...
    'Location','southoutside','NumColumns',numel(legendLabels),'FontSize',9);
set(leg,'Units','normalized','Position',[.08 .01 .84 .07]);
end
