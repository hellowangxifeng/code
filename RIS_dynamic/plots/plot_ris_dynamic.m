function files = plot_ris_dynamic(matFile,options)
%PLOT_RIS_DYNAMIC Draw every figure from a saved results file (no computation).
%   files = plot_ris_dynamic()                    results/ris_dynamic.mat
%   files = plot_ris_dynamic(matFile)
%   files = plot_ris_dynamic(matFile,options)     see rd_plot_options for fields
%   Only the light structure R written by rd_run is read, so figures can be
%   restyled at any time without touching the optics checkpoints.
%     Fig1 irradiance at the RX             Fig5 OP, height x wind grid
%     Fig2 ABER / OP, main condition        Fig6 complex field vs Monte Carlo
%     Fig3 ocean EGG x atmospheric GG       Fig7 required P_t vs height
%     Fig4 ABER, height x wind grid
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
files = [rd_fig_irradiance(R,o);rd_fig_main(R,o);rd_fig_turbulence(R,o); ...
    rd_fig_grid(R,o,'aber');rd_fig_grid(R,o,'op');rd_fig_validation(R,o);rd_fig_required(R,o)];
fprintf('Figures written to %s:\n',o.outputFolder);
fprintf('  %s\n',files);
end
