%RUN_WRIS_IRRADIANCE 修改下方参数区，点击运行即可重新计算并画图。
% +z -> right-hand rotation about fixed world x -> fixed world y.
% d0 = [sind(ay)*cosd(ax); -sind(ax); cosd(ay)*cosd(ax)].
wris_setup;
c=wris_config('production'); % 'smoke': 1e4 photons; 'production': 1e6

%% 仿真参数：入射角、风速、风向和采样时刻
c.comparison.rotationXYDeg=[0 30;0 60;0 80;45 45;-60 30;30 80];
windSpeeds=[3 5];                   % m/s; any number of wind speeds
windDirectionDeg=50;                % propagation direction: +Y toward +X
c.comparison.conditions=[windSpeeds(:),repmat(windDirectionDeg,numel(windSpeeds),1)];
c.comparison.site='center';         % center / crest / trough / rising / falling
c.ris.center=[0;0];                 % shared initial horizontal RIS centre (m)
c.sampling.comparisonTime=20;        % elapsed snapshot time (saved dt sample)
c.ris.duration=20; c.ris.dt=0.05; c.ris.startTime=25;
c.geometry.distance=10;             % initial TX -> lower-face centre distance (m)
c.output=fullfile(c.root,'results','incidence_xy');

%% 绘图参数：字号、分页和显示光强；修改这些参数不会使物理缓存失效
c.plot.rowsPerFigure=3;             % six angles -> two pages
c.plot.tickFontSize=16;             % 坐标轴及色条的数字字号
c.plot.labelFontSize=17;
c.plot.titleFontSize=12;
c.plot.figureTitleFontSize=16;
c.plot.resolution=300;
c.plot.quantity='noTurbulence';      % I/H; use 'irradiance' for EGG snapshot I
c.plot.visible='off';               % files always saved; 'on' also opens windows
recompute=true;                     % false: 只载入已保存数据重画；true: 计算当前物理参数

%% Calculation and rendering are independent
wris_results=calculateOrLoad(c,recompute);
wris_figure_files=wris_plot(wris_results);
disp(wris_figure_files);

function s=calculateOrLoad(c,recompute)
if recompute
    s=wris_run(c);
else
    saved=load(fullfile(c.output,'wris_results.mat'),'s');
    s=saved.s;                     % keeps the ORIGINAL physical configuration
    s.config.plot=c.plot;
    s.config.output=c.output;
end
end
