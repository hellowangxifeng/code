function rp_stage2_worker_refresh()
%RP_STAGE2_WORKER_REFRESH Reload changed receiver functions on idle workers.
clear rp_fading_response rp_link_metrics rp_stage2_wave rp_stage2_distance rp_stage2_select
rp_setup();
end
