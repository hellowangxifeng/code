%% RIS_dynamic one-click entry: floating RIS vs rough-sea direct link (MC)
% Press F5. Compute and plotting are decoupled: the computation writes
% results/ris_dynamic.mat (plus checkpoints), the plots only read that file.
% See README.md for the model, the parameters and the expected outputs.
clear; clc;
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'algorithm'),fullfile(here,'plots'));
rd_setup();

%% ---------------- choices ----------------
profile   = 'production';   % 'smoke' = minutes-long end-to-end check
doCompute = true;           % false: only redraw figures from the saved MAT
doPlot    = true;

%% ---------------- parameters (defaults in algorithm/rd_config.m) ----------------
c = rd_config(profile);
% c.scene.distancesM     = [7.5 10];        % receiver heights (m)
% c.sea.windsMps         = [2 3 4 5 6];     % sea states (calm U = 0 is always added)
% c.dynamic.frameStepS   = 0.5;             % optical frame step (s)
% c.ris.duration         = 30;              % motion record per wind (s)
% c.dynamic.photons      = 5e5;             % MC packets per link and frame
% c.dynamic.samples      = 2048;            % receiver jitter/pointing samples per frame
% c.study.fovHalfAngleRad = 0.130;          % FOV half-angle (rad)
% c.surface.crestBlocking = true;           % direct rays re-entering a crest are lost
% c.training.enabled     = true;            % sea-aware width selection
% c.incidence.enabled    = true;            % incidence-angle sweep (Fig7)
% c.incidence.anglesDeg  = [30 35 40 44 47 49 51 54];

plotOptions = struct();                      % see plots/rd_plot_options.m, e.g.
% plotOptions.curveWinds    = [0 3 5];
% plotOptions.fixedPowerDBm = 30;
% plotOptions.formats       = {'png','pdf','fig'};

%% ---------------- run ----------------
if doCompute
    R = rd_run(c);
end
if doPlot
    plot_ris_dynamic(fullfile(c.output,'ris_dynamic.mat'),plotOptions);
end
