function d=rp_stage2_apply_refinements(c,d)
%RP_STAGE2_APPLY_REFINEMENTS Restore resolved boundary candidates on clean runs.
files=dir(fullfile(c.output,'refined_candidates',sprintf('z%g_candidate*.mat',d.distanceM)));
for k=1:numel(files)
    token=regexp(files(k).name,'candidate(\d+)','tokens','once'); index=str2double(token{1});
    old=load(fullfile(files(k).folder,files(k).name),'b','opticalSignature');
    assert(strcmp(old.opticalSignature,rp_stage2_optical_signature(c)) && ...
        isequal(old.b.phase,d.bank{index}.phase),'rp:RefinedProvenance','Stale refined candidate.');
    b=old.b;
    b.channel=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,rp_link_samples(c));
    b.designMetrics=rp_link_metrics(c,b.channel,c.study.designPowerW); d.bank{index}=b;
end
d=rp_stage2_select(c,d);
end
