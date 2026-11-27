function b=rp_stage2_wave(c,g,mc,a,family,width,margin)
%RP_STAGE2_WAVE Restartable actual constrained phase and coherent propagation.
folder=fullfile(c.output,'wave_cache'); if ~isfolder(folder), mkdir(folder); end
stem=sprintf('a%05.2f_z%04.1f_%s_w%03d_m%+05.0f',g.incidenceDeg, ...
    g.rxPlaneZ-g.pRis(3)-c.ris.thickness,family,round(width*1000),margin*1e4);
file=fullfile(folder,[stem '.mat']); opticalSignature=rp_stage2_optical_signature(c);
if isfile(file)
    old=load(file,'b','signature','opticalSignature');
    assert(isfield(old,'opticalSignature') && strcmp(old.opticalSignature,opticalSignature), ...
        'rp:StaleStage2','Stale stage-two optical cache.');
    b=old.b;
    if ~strcmp(old.signature,c.study.signature)
        b.channel=rp_aperture_channel(c,g,b.wave,b.spot.centerXY,rp_link_samples(c));
        b.designMetrics=rp_link_metrics(c,b.channel,c.study.designPowerW);
        signature=c.study.signature; save(file,'b','signature','opticalSignature','-v7.3');
    end
    return
end
p=rp_finite_phase(c,g,a,family,width,margin);
r=rp_trace(c,g,a,p); win=rp_receiver_window(c,g,a,p,r);
p.centerXY=win.centerXY; p.receiverHalfWidthM=max(win.halfWidthM,[.32 .32]);
fprintf('%s full wave\n',stem); elapsed=tic;
w=rp_wave(c,g,a,mc,p,r); f=rp_spot_stats(w);
h=rp_aperture_channel(c,g,w,f.centerXY,rp_link_samples(c));
m=rp_link_metrics(c,h,c.study.designPowerW);
b=struct('family',family,'targetWidthM',width,'margin',margin,'phase',p, ...
    'wave',w,'spot',f,'channel',h,'designMetrics',m,'elapsedSeconds',toc(elapsed));
signature=c.study.signature; save(file,'b','signature','opticalSignature','-v7.3');
fprintf('DONE %s: w=%.3f cm, axis=%.3f, ABER=%.6g OP=%.6g\n', ...
    stem,100*f.equivalentRadiusM,f.axisRatio,m.aber,m.op);
end
