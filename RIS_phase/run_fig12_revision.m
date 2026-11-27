%RUN_FIG12_REVISION Recompute or redraw only the three-cm weak/weak comparison.
rp_setup();
storedStudy=load(fullfile('plots','stage2','study.mat'),'s');
[revision,files]=rp_fig12_run(storedStudy.s);
rp_fig12_report(revision);
