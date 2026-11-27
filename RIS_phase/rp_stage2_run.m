function s=rp_stage2_run(c)
%RP_STAGE2_RUN Reproducible second-stage entry; resume signed wave caches.
if nargin<1, c=rp_stage2_config(); end
rp_setup(); rp_study_check(c); if ~isfolder(c.output), mkdir(c.output); end
if ~isfolder(c.cache.folder), mkdir(c.cache.folder); end
c.study.signature=rp_stage2_signature(c); save(fullfile(c.output,'configuration.mat'),'c');
addpath(fullfile(c.root,'audit','stage2')); verify_joint_fading(c);
rp_stage2_scene(c,c.study.distancesM(1));
pool=gcp('nocreate'); if isempty(pool), pool=parpool('Processes',3); end
f=parfevalOnAll(pool,@rp_setup,0); fetchOutputs(f);
jobs=parallel.FevalFuture.empty(0,3);
for j=1:3, jobs(j)=parfeval(pool,@rp_stage2_distance_complete,1,c,c.study.distancesM(j)); end
s=struct('config',c,'distances',{cell(1,3)});
for j=1:3
    [idx,d]=fetchNext(jobs); s.distances{idx}=d;
    save(fullfile(c.output,'nominal_checkpoint.mat'),'s','-v7.3');
end
save(fullfile(c.output,'study.mat'),'s','-v7.3');
rp_stage2_finish(s);
end
