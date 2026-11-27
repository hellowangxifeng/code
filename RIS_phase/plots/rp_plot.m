function files=rp_plot(s)
%RP_PLOT Regenerate all finite-width study figures from saved data.
arguments
    s (1,1) struct
end
data=load(fullfile(s.config.output,'stage2.mat'),'extra');
files=rp_stage2_plot(s,data.extra);
end
