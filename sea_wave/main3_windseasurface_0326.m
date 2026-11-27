% Pure-sea companion to run_RIS_wind_comparison; no RIS dynamics.
% One pass, one output per case; the corrected spectrum is shared by both.
params = get_default_sea_params();
params.t0 = 27;
params.seed = 600;
exp_cases = [9,30;9,45;12,30];
output_dir = fullfile(fileparts(mfilename('fullpath')),'FigOutput_SeaSurface');
if ~isfolder(output_dir), mkdir(output_dir); end
sim_data = cell(1,size(exp_cases,1));
for i = 1:size(exp_cases,1)
    waveInputs = struct('g',params.g,'fMinHz',params.f_min,'fMaxHz',params.f_max, ...
        'frequencyBins',params.Fre_segment,'directionBins',params.Dir_segment, ...
        'windSpeed',exp_cases(i,1),'windAngleDeg',exp_cases(i,2),'seed',params.seed);
    wave = RIS_wave_components(waveInputs);
    field = struct('wave',wave,'parameters',struct('startTime',params.t0));
    [X,Y,Z] = RIS_sea_surface(field, ...
        linspace(params.x_range(1),params.x_range(2),params.grid_density), ...
        linspace(params.y_range(1),params.y_range(2),params.grid_density),0);
    sim_data{i} = struct('X',X,'Y',Y,'Z',Z,'wave',wave);
    fig = figure('Color','w','Visible','off','Position',[50,50,880,720]);
    ax = nexttile(tiledlayout(fig,1,1,'Padding','compact'));
    surf(ax,X,Y,Z,'EdgeColor','none');
    xlim(ax,params.x_range); ylim(ax,params.y_range);
    zlim(ax,[-0.5,0.5]); clim(ax,[-0.5,0.5]); zticks(ax,[-0.5,0,0.5]);
    pbaspect(ax,[1,1,0.7]); view(ax,-37.5,30); grid(ax,'on'); box(ax,'on');
    colormap(ax,turbo(256)); cb = colorbar(ax); cb.Label.String = 'Elevation (m)';
    xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); zlabel(ax,'z (m)');
    title(ax,sprintf('U=%g m/s | wind %g deg | t_{sea}=%g s',exp_cases(i,:),params.t0));
    name = sprintf('Surf_U%02d_Phi%02d',exp_cases(i,:));
    savefig(fig,fullfile(output_dir,[name,'.fig']));
    exportgraphics(fig,fullfile(output_dir,[name,'.png']),'Resolution',160);
    close(fig);
end
save(fullfile(output_dir,'sea_surface_comparison.mat'),'sim_data','params','exp_cases');
