%% RIS_dynamic one-click entry: frozen rough sea, floating RIS vs no-RIS link
% Press F5. Compute and plotting are decoupled: the computation writes
% results/ris_dynamic.mat (plus checkpoints and CSV tables), the plots only
% read that file. See 算法与图件说明.md for the model and every figure.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'algorithm'),fullfile(here,'plots'));
rd_setup();

%% ---------------- choices ----------------
profile   = 'production';   % 'smoke' = short end-to-end check (results/smoke)
doCompute = true;           % false: only redraw figures from the saved MAT
doPlot    = true;

%% ---------------- parameters (defaults in algorithm/rd_config.m) ----------------
c = rd_config(profile);
% c.scene.distancesM     = [10 20 30 50];  % RX plane above the calm RIS centre (m)
% c.sea.windsMps         = [3 9 12];       % sea states, 3-25 m/s
% c.main.distanceM       = 30;             % Fig1-3, Fig6 condition (must be in the grids)
% c.main.windMps         = 9;
% c.sea.freezeTimeS      = 15;             % frozen sea instant (s)
% c.design.widthsM       = [.04 .06 .08 .10 .13 .16 .20 .25 .30];  % Flat/Robust target widths
% c.study.fovHalfAngleRad = 0.130;         % FOV half-angle (rad)
% c.validation.photons   = 1e6;            % MC check of the complex field
% c.direct.photons       = 2e6;            % no-RIS rough-surface MC
% c.metrics.policy       = 'per_power';    % or 'design_power' (one phase for all powers)
% c.run.parallel         = true;           % parfor over the designs of one condition

plotOptions = struct();                      % see plots/rd_plot_options.m, e.g.
% plotOptions.formats       = {'png','pdf','fig'};
% plotOptions.mapHalfWidthM = 0.5;

%% ---------------- run ----------------
if doCompute
    R = rd_run(c);
end
if doPlot
    plot_ris_dynamic(fullfile(c.output,'ris_dynamic.mat'),plotOptions);
end
