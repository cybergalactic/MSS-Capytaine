function geometry = plotVesselGeometry(case_name)
% plotVesselGeometry compares source geometry and saved Capytaine hull panels.
% Top: CSV/MGF sections or GDF/STL source panels. Bottom: saved solver mesh.
% Before solving, the bottom subplot explains that no saved mesh is available.
% Both use MSS FSD; no Python process or MSS installation is needed to plot.
% Examples:
%   plotVesselGeometry('testShip_csv')
%   plotVesselGeometry('otter_80kg_csv')
%   plotVesselGeometry('supply_mgf')
%   plotVesselGeometry('tanker_gdf')
%   plotVesselGeometry('wigley_stl')
% A JSON file path is also accepted. If the geometry filename has no known
% extension, the reader identifies its format from the file contents.
% MGF uses draft_m from JSON. STL must already be in metres and MSS FSD.
% CSV demihull_spacing_m places two identical hulls at y = +/-spacing/2.
% The mesh is read from output_dir/generated_hull_panels.mat (default results/).
% Plotting never generates panels or reruns the solver. Rerun Python after
% changing geometry or meshing settings to refresh the saved mesh.
geometry=readVesselGeometry(case_name);
geometry.mesh=readCapytaineMesh(geometry);
tag=['MSSCapytaineGeometry:' geometry.path];
fig=findobj(0,'Type','figure','Tag',tag);
if isempty(fig)
    fig=figure('Name',[geometry.name ' geometry'],'NumberTitle','off','Tag',tag);
else
    fig=fig(1); clf(fig); figure(fig);
end
ax=subplot(2,1,1,'Parent',fig); hold(ax,'on');
if isempty(geometry.faces)
    if strcmp(geometry.format,'mgf')
        drawSections(ax,geometry.sections,[0.75 0.75 0.75],0.5,geometry.section_centers_m);
        drawSections(ax,geometry.wetted_sections,[0.05 0.35 0.7],1.0,geometry.section_centers_m);
        subtitle='MGF sections: blue wetted, grey above water';
    else
        drawSections(ax,geometry.sections,[0.05 0.35 0.7],0.8,geometry.section_centers_m);
        subtitle='CSV offset sections';
        if numel(geometry.section_centers_m)==2
            subtitle='CSV demi-hull sections, placed at configured spacing';
        end
    end
else
    drawPanels(ax,geometry.vertices,geometry.faces);
    subtitle=[upper(geometry.format) ' source panels'];
end
drawCG(ax,geometry.center_of_mass_m);
finishAxes(ax,{[geometry.name ' source geometry'],subtitle});
mesh_ax=subplot(2,1,2,'Parent',fig);
if geometry.mesh.available
    hold(mesh_ax,'on');
    drawPanels(mesh_ax,geometry.mesh.vertices,geometry.mesh.faces);
    drawCG(mesh_ax,geometry.mesh.center_of_mass_m);
    finishAxes(mesh_ax,{'Saved Capytaine hull mesh', ...
        sprintf('%d panels (MSS FSD)',size(geometry.mesh.faces,1))});
    % Keep both hulls at the same physical scale, including dry MGF sections.
    for property={'XLim','YLim','ZLim'}
        limits=[get(ax,property{1});get(mesh_ax,property{1})];
        set([ax mesh_ax],property{1},[min(limits(:,1)),max(limits(:,2))]);
    end
    % Retain the link so rotations and axis-limit changes stay synchronized.
    axes_link=linkprop([ax mesh_ax],{'View','XLim','YLim','ZLim'});
    setappdata(fig,'GeometryAxesLink',axes_link);
else
    axis(mesh_ax,'off');
    title(mesh_ax,'Capytaine hull mesh');
    text(mesh_ax,0.5,0.5,{'No saved mesh yet.', ...
        'Run Python for this case, then plot again.'}, ...
        'Units','normalized','HorizontalAlignment','center','Interpreter','none');
end
end

function drawPanels(ax,vertices,faces)
% NaN padding lets patch preserve quadrilaterals and repeated-vertex triangles.
faces=double(faces);
if size(faces,2)==4
    triangles=faces(:,3)==faces(:,4); faces(triangles,4)=NaN;
end
patch('Parent',ax,'Faces',faces,'Vertices',vertices, ...
    'FaceColor',[0.65 0.78 0.9],'EdgeColor',[0.2 0.3 0.4],'LineWidth',0.25);
end

function drawCG(ax,cg)
if isempty(cg), return; end
plot3(ax,cg(1),cg(2),cg(3),'rx','MarkerSize',9,'LineWidth',1.5);
text(ax,cg(1),cg(2),cg(3),' CG','Color',[0.8 0 0]);
end

function finishAxes(ax,heading)
hold(ax,'off'); axis(ax,'equal'); grid(ax,'on'); view(ax,3);
xlabel(ax,'x forward [m]'); ylabel(ax,'y starboard [m]'); zlabel(ax,'z down [m]');
set(ax,'ZDir','reverse'); title(ax,heading,'Interpreter','none');
end

function drawSections(ax,sections,color,width,centers)
for center=centers
    for i=1:numel(sections)
        points=sections{i};
        plot3(ax,points(:,1),center+points(:,2),points(:,3),'-','Color',color,'LineWidth',width);
        plot3(ax,points(:,1),center-points(:,2),points(:,3),'-','Color',color,'LineWidth',width);
    end
end
end
