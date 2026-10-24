function exportGeometryPreview()
% Export README artwork from the public MATLAB geometry viewer; no BEM solve.
% From the repository root: addpath('tests'); exportGeometryPreview()
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'matlab'));
output = fullfile(root,'docs','images');
if ~isfolder(output), mkdir(output); end
visibility = get(groot,'DefaultFigureVisible');
set(groot,'DefaultFigureVisible','off');
restore = onCleanup(@() set(groot,'DefaultFigureVisible',visibility)); %#ok<NASGU>
overview = figure('Visible','off','Color','white', ...
    'Units','pixels','Position',[50 50 2000 1100]);
cleanup = onCleanup(@() close(overview)); %#ok<NASGU>
cases = {'testShip','LAUV_marie'};
labels = {'Surface vessel: testShip','Submerged vehicle: LAUV Marie'};
for column = 1:numel(cases)
    plotVesselGeometry(cases{column});
    source = gcf;
    axesList = findall(source,'Type','axes');
    heights = arrayfun(@(ax) ax.Position(2),axesList);
    [~,order] = sort(heights,'descend');
    for row = 1:2
        ax = copyobj(axesList(order(row)),overview);
        ax.Units = 'normalized';
        % Reserve room outside the axes for the long FSD labels in 3-D.
        ax.Position = [0.10+(column-1)*0.5, 0.57-(row-1)*0.47, 0.32, 0.33];
        ax.FontName = 'Arial';
        ax.FontSize = 12;
        ax.GridAlpha = 0.15;
        ax.LineWidth = 0.75;
        if row == 1
            title(ax,{labels{column},'CSV offset sections'},'FontSize',15,'FontWeight','normal');
        else
            title(ax,'Saved Capytaine hull mesh','FontSize',15,'FontWeight','normal');
        end
        ax.Title.Interpreter = 'none';
    end
    close(source);
end
exportgraphics(overview,fullfile(output,'vessel-geometry.png'), ...
    'Resolution',150,'BackgroundColor','white','Padding','figure');
fprintf('Exported docs/images/vessel-geometry.png from plotVesselGeometry.\n');
end
