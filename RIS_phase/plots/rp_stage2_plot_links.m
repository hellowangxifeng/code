function files=rp_stage2_plot_links(d,c)
%RP_STAGE2_PLOT_LINKS Common fixed receiver; per-power optimal loaded phases.
data=rp_fixed_receiver_links(d,c,'per_power');
save(fullfile(c.output,'communication_data.mat'),'data');
writetable(data.curveTable,fullfile(c.output,'communication_curves.csv'));
files={rp_stage2_plot_link_metric(data,c,'aber'),rp_stage2_plot_link_metric(data,c,'op')};
end
