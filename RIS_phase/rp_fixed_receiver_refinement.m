function refined=rp_fixed_receiver_refinement(c,d,b)
%RP_FIXED_RECEIVER_REFINEMENT Resolve any selected core near the 2 cm bound.
rp_setup();file=fullfile(c.output,'fixed_receiver',sprintf('refined_%s_w%03d.mat',b.family,round(1000*b.targetWidthM)));
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native([rp_stage2_optical_signature(c),jsonencode(b.phase), ...
    jsonencode(c.study),fileread([mfilename('fullpath') '.m'])],'UTF-8'));
signature=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
if isfile(file)
    old=load(file,'refined','signature');
    if strcmp(old.signature,signature),refined=old.refined;return;end
end
cc=c;cc.spot.mapBins=[192 192];cc.wave.apertureOversampling=2;cc.wave.fftPadding=4;
r=rp_trace(cc,d.geometry,d.illumination,b.phase);
w=rp_wave(cc,d.geometry,d.illumination,d.mc,b.phase,r); f=rp_spot_stats(w);
samples=rp_link_samples(c,c.study.validationSamples,c.study.seed+703);
baseH=rp_aperture_channel(c,d.geometry,b.wave,b.receiverCenterXY,samples);
newH=rp_aperture_channel(c,d.geometry,w,b.receiverCenterXY,samples);
baseM=rp_fixed_receiver_metrics(c,baseH,1);newM=rp_fixed_receiver_metrics(c,newH,1);
check=struct('widthRelativeChange',abs(f.equivalentRadiusM/b.spot.equivalentRadiusM-1), ...
    'powerRelativeChange',abs(newM.meanReceivedPowerW/baseM.meanReceivedPowerW-1), ...
    'aberAbsoluteChange',abs(newM.aber-baseM.aber),'opAbsoluteChange',abs(newM.op-baseM.op), ...
    'refinedR50M',f.r50M);
assert(check.widthRelativeChange<.02 && check.powerRelativeChange<.03 && ...
    check.aberAbsoluteChange<.002 && check.opAbsoluteChange<.005,'rp:FixedCoreConvergence','Fixed-RX core needs more resolution.');
b.wave=w;b.spot=f;b.channel=rp_aperture_channel(c,d.geometry,w,b.receiverCenterXY,rp_link_samples(c));
b.designMetrics=rp_fixed_receiver_metrics(c,b.channel,c.study.designPowerW);b.refinement=check;
refined=struct('b',b,'check',check);save(file,'refined','signature','-v7.3');
end
