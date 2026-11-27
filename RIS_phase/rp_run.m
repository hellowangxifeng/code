function s=rp_run(c)
%RP_RUN Public entry for the finite-width, three-distance receiver study.
arguments
    c (1,1) struct = rp_stage2_config()
end
s=rp_stage2_run(c);
end
