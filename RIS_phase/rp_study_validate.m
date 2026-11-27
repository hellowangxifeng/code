function v=rp_study_validate(s)
%RP_STUDY_VALIDATE Compatibility wrapper for the current stage-two audit.
data=load(fullfile(s.config.output,'stage2.mat'),'extra');
v=rp_stage2_validate(s,data.extra);
end
