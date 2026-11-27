function files = plot_wris_dynamic(resultFile,opts)
%PLOT_WRIS_DYNAMIC Draw every figure from SAVED results only.
%   files = plot_wris_dynamic()                 reads ../results/wris_dynamic.mat
%   files = plot_wris_dynamic(resultFile)       any wris_dynamic.mat (e.g. smoke)
%   files = plot_wris_dynamic(resultFile,opts)  opts fields (all optional):
%       visible      'off' (default) | 'on'
%       resolution   PNG dpi, default 300
%       fontName     default 'Times New Roman'
%       fontSize     default 10
%       mapTimeS     frame time shown in the irradiance maps, default 10 s
%       mapCondition wind-condition column for the maps, default = last
%       outputFolder default = folder of resultFile
%   No MC transport, motion solve or algorithm function is called here, so
%   fonts/layout can be changed and re-run in seconds.
%
%   Figures (PNG + editable FIG):
%     fig1_rx_power_time   RX power vs time, all links, same sea and RX
%     fig2_rx_power_cdf    empirical CDF of RX power over frames + outage threshold
%     fig3_sea_tir         rough vs still-water TIR share and sea slope at the board
%     fig4_rx_maps         detector-plane irradiance at one common time
%     fig5_summary_wind    mean/min/max RX power and outage vs wind speed
here = fileparts(mfilename('fullpath'));
if nargin<1 || isempty(resultFile)
    resultFile = fullfile(fileparts(here),'results','wris_dynamic.mat');
end
if nargin<2, opts = struct(); end
resultFile = char(resultFile);
assert(isfile(resultFile),'wrd:NoResults', ...
    'Result file not found: %s\nRun run_wris_dynamic first (or main_ris_dynamic).',resultFile);
addpath(here);
loaded = load(resultFile,'results');
r = loaded.results;
r.config = defaults(r.config);
opts = options(opts,resultFile,size(r.cases,2));
if ~isfolder(opts.outputFolder), mkdir(opts.outputFolder); end
files = [wrd_fig_power(r,opts); wrd_fig_cdf(r,opts); wrd_fig_tir(r,opts); ...
    wrd_fig_maps(r,opts); wrd_fig_summary(r,opts)];
fprintf('Figures written to %s\n',opts.outputFolder);
fprintf('  %s\n',files(endsWith(files,'.png')));
end

function opts = options(opts,resultFile,nConditions)
base = struct('visible','off','resolution',300,'fontName','Times New Roman', ...
    'fontSize',10,'mapTimeS',10,'mapCondition',nConditions,'outputFolder',fileparts(resultFile));
names = fieldnames(opts);
for j = 1:numel(names)
    assert(isfield(base,names{j}),'wrd:PlotOption','Unknown plot option: %s',names{j});
    base.(names{j}) = opts.(names{j});
end
opts = base;
opts.visible = validatestring(opts.visible,{'on','off'});
if isempty(opts.outputFolder), opts.outputFolder = fileparts(resultFile); end
if isempty(opts.outputFolder), opts.outputFolder = pwd; end
opts.mapCondition = min(max(1,round(opts.mapCondition)),nConditions);
end

function c = defaults(c)
% Results written by the first RIS_dynamic version lack these fields.
if ~isfield(c,'metrics'), c.metrics = struct(); end
if ~isfield(c.metrics,'outageThresholdW'), c.metrics.outageThresholdW = 1e-6; end
if ~isfield(c.metrics,'displayFloorW'), c.metrics.displayFloorW = 1e-9; end
end
