function d=rp_stage2_promote(c,d,index,check)
%RP_STAGE2_PROMOTE Use resolved optical data when a boundary core fails.
b=d.bank{index}; b.wave=check.refinedWave; b.spot=check.refinedSpot;
b.channel=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,rp_link_samples(c));
b.designMetrics=rp_link_metrics(c,b.channel,c.study.designPowerW);
opticalSignature=rp_stage2_optical_signature(c);
folder=fullfile(c.output,'refined_candidates'); if ~isfolder(folder), mkdir(folder); end
save(fullfile(folder,sprintf('z%g_candidate%02d.mat',d.distanceM,index)), ...
    'b','opticalSignature','-v7.3');
d.bank{index}=b; d=rp_stage2_select(c,d);
end
