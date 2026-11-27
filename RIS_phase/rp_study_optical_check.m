function check=rp_study_optical_check(c,d,index)
%RP_STUDY_OPTICAL_CHECK Refined optical cache; refresh current EGG/GG statistics.
if nargin<3, index=d.objectiveIndices(1); end
file=fullfile(c.output,sprintf('optical_check_z%g_candidate%02d.mat',d.distanceM,index));
item=d.bank{index}; digest=java.security.MessageDigest.getInstance('SHA-256');
key=[rp_stage2_optical_signature(c),jsonencode(item.phase), ...
    jsonencode([c.study.waveCheckBins 2 4]),fileread([mfilename('fullpath') '.m'])];
digest.update(unicode2native(key,'UTF-8'));
signature=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
reuse=false;
if isfile(file)
    data=load(file,'check','signature');
    if strcmp(data.signature,signature), w=data.check.refinedWave; reuse=true; end
end
if ~reuse
    refined=c; refined.spot.mapBins=[c.study.waveCheckBins c.study.waveCheckBins];
    refined.wave.apertureOversampling=2; refined.wave.fftPadding=4;
    r=rp_trace(refined,d.geometry,d.illumination,item.phase);
    w=rp_wave(refined,d.geometry,d.illumination,d.mc,item.phase,r);
end
f=rp_spot_stats(w); samples=rp_link_samples(c,c.study.validationSamples,c.study.seed+703);
oldH=rp_aperture_channel(c,d.geometry,item.wave,item.spot.centerXY,samples);
newH=rp_aperture_channel(c,d.geometry,w,f.centerXY,samples);
oldM=rp_link_metrics(c,oldH,c.study.designPowerW,64);
newM=rp_link_metrics(c,newH,c.study.designPowerW,64);
check=struct('distanceM',d.distanceM,'candidate',index,'refinedWave',w, ...
    'refinedSpot',f,'refinedMetrics',newM,'baseMetrics',oldM, ...
    'relativeWidthChange',abs(f.equivalentRadiusM/item.spot.equivalentRadiusM-1), ...
    'relativePowerChange',abs(newM.meanReceivedPowerW/oldM.meanReceivedPowerW-1), ...
    'centroidChangeM',norm(f.centerXY-item.spot.centerXY), ...
    'aberAbsoluteChange',abs(newM.aber-oldM.aber),'opAbsoluteChange',abs(newM.op-oldM.op), ...
    'opticalSignature',signature,'linkSignature',c.study.signature);
check.passed=check.relativeWidthChange<.02 && check.relativePowerChange<.03 && ...
    check.centroidChangeM<.002 && check.aberAbsoluteChange<.002 && check.opAbsoluteChange<.005;
save(file,'check','signature','-v7.3');
assert(check.passed,'rp:OpticalConvergence','Optical refinement failed at z=%g m.',d.distanceM);
end
