function files = plot_ris_dynamic(matFile,options)
%PLOT_RIS_DYNAMIC Draw every figure from a saved results file (no computation).
%   files = plot_ris_dynamic()                    results/ris_dynamic.mat
%   files = plot_ris_dynamic(matFile)
%   files = plot_ris_dynamic(matFile,options)     see rd_plot_options for fields
%   Only the light structure R written by rd_run is read, so figures can be
%   restyled at any time without touching the Monte Carlo checkpoints.
here = fileparts(mfilename('fullpath'));
if nargin<1 || isempty(matFile)
    matFile = fullfile(fileparts(here),'results','ris_dynamic.mat');
end
if nargin<2, options = struct(); end
assert(isfile(matFile),'rd:Plot','Results file not found: %s (run rd_run first).',matFile);
loaded = load(matFile,'R');
R = loaded.R;
o = rd_plot_options(R,options,matFile);
if ~isfolder(o.outputFolder), mkdir(o.outputFolder); end
files = strings(0,1);
files = [files;rd_fig_irradiance(R,o)];
files = [files;rd_fig_curves(R,o,'aber')];
files = [files;rd_fig_curves(R,o,'op')];
files = [files;rd_fig_timeseries(R,o)];
files = [files;rd_fig_cdf(R,o)];
files = [files;rd_fig_wind(R,o)];
if isstruct(R.incidence) && ~isempty(R.incidence)
    files = [files;rd_fig_incidence(R,o)];
end
if height(R.trainingTable)>0
    files = [files;rd_fig_width(R,o)];
end
fprintf('Figures written to %s:\n',o.outputFolder);
fprintf('  %s\n',files);
end
