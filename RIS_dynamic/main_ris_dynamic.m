%MAIN_RIS_DYNAMIC 一键入口：动态海况下 浮动RIS 与 粗糙海面直达链路(MC) 对比。
%   用法：在 MATLAB 中打开本文件，按需修改“参数区”，点击“运行”(F5)。
%   1) 计算：run_wris_dynamic(c)  -> results/wris_dynamic.mat + CSV + run_log.txt
%   2) 画图：plot_wris_dynamic(...) -> results/fig*.png + fig*.fig（只读 MAT）
%   计算与画图完全解耦：改字号/画幅只需 doCompute=false 再运行一次。
%   详细说明见同目录 README.md。

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'algorithm'),fullfile(here,'plots'));
wrd_setup();

%% ===================== 参数区 =====================
profile   = 'production';   % 'smoke'：1个风况×1个角度×1 s×2000光子，先用它检查环境
doCompute = true;           % false：跳过计算，只用已保存的 MAT 重画
doPlot    = true;

c = wrd_config(profile);    % 全部默认参数见 algorithm/wrd_config.m
% ---- 以下为可选覆盖示例（取消注释即可；smoke 档建议不要改）----
% c.comparison.conditions    = [3 50; 5 50];   % [风速 m/s, 风向 deg]，同一海况驱动 RIS 与直达链路
% c.comparison.rotationXYDeg = [0 40; 0 49];   % TX 光轴：+z 先绕世界x转、再绕世界y转 (deg)
% c.comparison.frameStep     = 0.5;            % 输出帧间隔 s（必须是 c.ris.dt=0.05 s 的整数倍）
% c.comparison.flatBaseline  = true;           % 同时计算旧的“静止海面直达”基线
% c.comparison.risWaveBlocking = false;        % true：RIS->RX 光线也考虑被波峰遮挡（敏感性分析）
% c.mc.photons               = 100000;         % 每条链路每帧光子包数
% c.receiver.center          = [8.3 0 5];      % 固定接收机位置 (m)
% c.receiver.radius          = 0.05;           % 接收孔径半径 (m)
% c.receiver.fovDeg          = 15;             % 视场半角 (deg)
% c.metrics.outageThresholdW = 1e-6;           % 中断门限 (W)，按实际接收机灵敏度设置
% c.run.resume               = true;           % 中断后重跑会复用已完成风况的 checkpoint

plotOptions = struct('visible','off', ...       % 'on' 同时弹出窗口
    'resolution',300,'fontName','Times New Roman','fontSize',10, ...
    'mapTimeS',10);                             % 光强分布图取的时刻 (s)
%% ==================================================

resultFile = fullfile(c.output,'wris_dynamic.mat');
if doCompute
    results = run_wris_dynamic(c); %#ok<NASGU>
end
if doPlot
    figureFiles = plot_wris_dynamic(resultFile,plotOptions);
end
