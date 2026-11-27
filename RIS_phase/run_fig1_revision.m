%RUN_FIG1_REVISION Rebuild Fig1 from the saved common lower-surface MC packets.
rp_setup();root=fileparts(mfilename('fullpath'));addpath(fullfile(root,'audit','fig1_revision'));
saved=load(fullfile(root,'plots','stage2','study.mat'),'s');
rp_stage2_plot_irradiance(saved.s);
fig1Verification=verify_fig1_revision(saved.s);
disp(fig1Verification.status);
