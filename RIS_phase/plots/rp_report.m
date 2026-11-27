function rp_report(s)
%RP_REPORT Export finite-width receiver tables and actual phase commands.
arguments
    s (1,1) struct
end
data=load(fullfile(s.config.output,'stage2.mat'),'extra');
validation=load(fullfile(s.config.output,'validation.mat'),'v');
rp_stage2_report(s,data.extra,validation.v);
end
