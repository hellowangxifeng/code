function [s,extra,validation]=rp_stage2_finish(s)
%RP_STAGE2_FINISH All additional second-stage scans, validation and figures.
c=s.config; pool=gcp('nocreate'); if isempty(pool), pool=parpool('Processes',3); end
% Resolve a narrow low-jitter core before accepting its feasibility.
while true
    extra=rp_stage2_receiver(s); low=extra.widthJitter.selection;
    id=low.candidate(find(low.objectiveId==1 & low.sigmaM==0,1));
    core=rp_study_optical_check(c,s.distances{2},id);
    if core.refinedSpot.r50M>=.02, break; end
    s.distances{2}=rp_stage2_promote(c,s.distances{2},id,core);
end
save(fullfile(c.output,'study.mat'),'s','-v7.3');
angleJobs=parallel.FevalFuture.empty(0,2*numel(c.study.angleGridDeg)); count=0;
for j=1:2
    design=rp_stage2_compact_design(s.distances{j+1});
    for k=1:numel(c.study.angleGridDeg)
        count=count+1; angleJobs(count)=parfeval(pool,@rp_stage2_angle_case,1,c,design,k);
    end
end
opticalJobs=parallel.FevalFuture.empty(0,6);
for j=1:3, opticalJobs(j)=parfeval(pool,@rp_study_optical_check,1,c,s.distances{j}); end
opticalJobs(4)=parfeval(pool,@rp_study_optical_check,1,c,s.distances{3},1);
opticalJobs(5)=parfeval(pool,@rp_study_optical_check,1,c,s.distances{1},s.distances{1}.objectiveIndices(3));
opticalJobs(6)=parfeval(pool,@rp_study_optical_check,1,c,s.distances{2},id);
extra.angles=cell(1,2); extra.optical=cell(1,6);
rp_study_progress(c,'incidence_and_optical_refinement',0,numel(angleJobs)+numel(opticalJobs));
cases=cell(1,numel(angleJobs));
for j=1:numel(angleJobs)
    [idx,value]=fetchNext(angleJobs); cases{idx}=value;
    rp_study_progress(c,'incidence_and_optical_refinement',j,numel(angleJobs)+numel(opticalJobs));
end
for j=1:2, extra.angles{j}=rp_stage2_merge_angles(c,cases,c.study.distancesM(j+1)); end
for j=1:numel(opticalJobs)
    [idx,value]=fetchNext(opticalJobs); extra.optical{idx}=value;
    rp_study_progress(c,'incidence_and_optical_refinement',numel(angleJobs)+j,numel(angleJobs)+numel(opticalJobs));
end
save(fullfile(c.output,'stage2.mat'),'s','extra','-v7.3');
rp_study_progress(c,'independent_validation'); validation=rp_stage2_validate(s,extra);
rp_study_progress(c,'plotting_and_export');
rp_stage2_plot(s,extra); graphicsAudit=rp_stage2_graphics_audit(c);
rp_stage2_report(s,extra,validation);
save(fullfile(c.output,'stage2.mat'),'s','extra','validation','graphicsAudit','-v7.3');
rp_study_progress(c,'numerical_complete_visual_review_pending',1,1,'review');
end
