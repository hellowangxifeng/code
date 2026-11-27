function RIS_draw_state(ax,motion,time,scene,heightLimits)
%RIS_DRAW_STATE Draw the configured glass plane on its incident sea snapshot.
% The load-bearing frame participates in dynamics but is not drawn.
% scene has X,Y,Z in global metres; time is an existing elapsed sample.
    arguments
        ax (1,1) matlab.graphics.axis.Axes
        motion (1,1) struct
        time (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
        scene (1,1) struct
        heightLimits (1,2) double {mustBeReal,mustBeFinite} = [-0.5,0.5]
    end
    index = find(abs(motion.t-time)<1e-10,1);
    if isempty(index)
        error('RIS_draw_state:Time','Choose an existing sample of motion.t.');
    end
    surf(ax,scene.X,scene.Y,scene.Z,'EdgeColor','none');
    hold(ax,'on');
    R = motion.rotationMatrix(:,:,index);
    patch(ax,'Vertices',motion.corners(:,:,index),'Faces',[1,2,3,4], ...
        'FaceColor',[0.45,0.85,0.95],'FaceAlpha',1, ...
        'EdgeColor','k','LineWidth',1.3,'Tag','RISGlassPlane');
    normal = motion.normal(:,index); center = motion.boardCenter(:,index);
    quiver3(ax,center(1),center(2),center(3),0.5*normal(1),0.5*normal(2),0.5*normal(3), ...
        0,'k','LineWidth',1.4);
    colors = [0.8,0.12,0.12;0.1,0.55,0.2];
    for bodyAxis = 1:2
        direction = 0.4*min(motion.parameters.Lx,motion.parameters.Ly)*R(:,bodyAxis);
        quiver3(ax,center(1),center(2),center(3),direction(1),direction(2),direction(3), ...
            0,'Color',colors(bodyAxis,:),'LineWidth',1.5);
    end
    xlim(ax,[min(scene.X(:)),max(scene.X(:))]);
    ylim(ax,[min(scene.Y(:)),max(scene.Y(:))]);
    zlim(ax,heightLimits); clim(ax,heightLimits);
    % Fixed display proportions; numerical heights are never rescaled.
    pbaspect(ax,[1,1,0.7]); view(ax,-37.5,30); grid(ax,'on'); box(ax,'on');
    zticks(ax,[-0.5,0,0.5]);
    colormap(ax,turbo(256)); cb = colorbar(ax); cb.Label.String = 'Elevation (m)';
    xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); zlabel(ax,'z (m)'); set(ax,'FontSize',9);
end
