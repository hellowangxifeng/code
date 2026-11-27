function v=rp_validate(s)
%RP_VALIDATE Validate finite-width candidate power, phase and receiver statistics.
arguments
    s (1,1) struct
end
data=load(fullfile(s.config.output,'stage2.mat'),'extra');
v=rp_stage2_validate(s,data.extra);
end
