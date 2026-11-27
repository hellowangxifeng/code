function e = rd_candidate(c,s,g,target,spec,samples,rx)
%RD_CANDIDATE One RIS design on the frozen sea: complex field AND full MC.
%   e = rd_candidate(c,s,g,target,spec,samples,rx)
%   Complex field (the formula used for all curves):
%     phase p (rd_phase) -> scattered MC packets through p (rp_trace) ->
%     receiver window (rp_receiver_window, widened to cover the RX jitter)
%     -> ballistic Fresnel/CZT field + scattered halo (rd_wave) ->
%     tilted-disc aperture integral per receiver sample (rp_aperture_channel).
%   Monte Carlo (independent check): an independent wris_mc run of the same
%   frozen lower face (s.lcp, raw packet powers, ballistic AND scattered),
%   traced through the same loaded gradients (rd_trace), captured by the
%   same physical disc samples (rd_capture, per-ray FOV) and binned on the
%   same pixels (rd_map).
p = rd_phase(c,g,s.a,s.base,spec.family,spec.width,target);
r = rp_trace(c,g,s.a,p);
native = rp_receiver_window(c,g,s.a,p,r);
half = max(c.wave.receiverHalfWidthM,abs(native.centerXY-target)+native.halfWidthM);
pixel = c.rx.apertureRadius/c.wave.pixelsPerRadius;
bins = min(c.wave.maxBins,max(c.wave.minBins,ceil(2*half/pixel)));
cw = c;
cw.spot.mapBins = bins;
p.centerXY = target;
p.receiverHalfWidthM = half;
clock = tic;
w = rd_wave(cw,g,s.a,s.mc,p,r);
spot = rp_spot_stats(w);
h = rp_aperture_channel(c,g,w,target,samples);

land = rd_trace(c,g,s.lcp,p,g.rxPlaneZ);
gainMC = rd_capture(c,rd_compress(land,c.compress),rx);
dx = w.x(2)-w.x(1);
dy = w.y(2)-w.y(1);
mapMC = rd_map(land,[w.x-dx/2,w.x(end)+dx/2],[w.y-dy/2,w.y(end)+dy/2]);

limit = floor(c.device.maxAngleDeg/c.device.angleStepDeg)*c.device.angleStepDeg;
legal = max(abs(p.commandDeg),[],'all')<=limit+1e-12 && ...
    max(abs(p.commandDeg/c.device.angleStepDeg-round(p.commandDeg/c.device.angleStepDeg)),[],'all')<1e-9;
assert(legal,'rd:Phase','Illegal LC command in %s design.',spec.family);
e = struct('family',spec.family,'width',spec.width,'targetXY',target, ...
    'gainField',single(h.withoutFOV),'gainMC',single(gainMC), ...
    'outsideMapFraction',h.outsideMapApertureFraction, ...
    'x',w.x,'y',w.y,'mapField',single(w.irradianceWm2),'mapMC',single(mapMC), ...
    'mapScatter',single(w.scatteredIrradianceWm2),'spot',spot, ...
    'maxCommandDeg',max(abs(p.commandDeg),[],'all'),'saturatedPowerFraction',p.saturatedPowerFraction, ...
    'tiles',numel(p.unitIds),'pupilPowerW',w.pupilPowerW,'mapPowerW',w.mapPowerW, ...
    'scatteredMapPowerW',w.scatteredMapPowerW,'beamletTruncatedPowerW',w.beamletTruncatedPowerW, ...
    'omittedBallisticPowerW',w.omittedBallisticPowerW,'coherenceRingPowerW',w.coherenceRingPowerW, ...
    'mcPlanePowerW',sum(land.power),'mcMapPowerW',sum(mapMC,'all')*dx*dy, ...
    'mcMapBallisticPowerW',sum(land.power(land.ballistic & ...
        abs(land.xy(:,1)-target(1))<=half(1) & abs(land.xy(:,2)-target(2))<=half(2))), ...
    'elapsedSeconds',toc(clock));
fprintf('    %-8s w = %5.3f m | R_eq %6.2f cm | map %.4g W (MC %.4g W) | <G> field %.4g, MC %.4g | %.0f s\n', ...
    spec.family,spec.width,100*spot.equivalentRadiusM,w.mapPowerW,e.mcMapPowerW, ...
    mean(h.withoutFOV),mean(gainMC),e.elapsedSeconds);
end
