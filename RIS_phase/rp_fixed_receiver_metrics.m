function m=rp_fixed_receiver_metrics(c,h,powerW)
%RP_FIXED_RECEIVER_METRICS Exact Rice tail with conditional aperture statistics.
accepted=h.angleErrorRad<=c.study.fovHalfAngleRad & h.angleErrorRad<pi/2;
conditional=struct('gain',h.withoutFOV(accepted),'fovOutageExact',0,'fovOutageSample',0);
assert(any(accepted),'rp:FixedFOV','No in-FOV samples.');
m=rp_link_metrics(c,conditional,powerW,64);
if c.study.angleSigmaRad>0
    pOut=ncx2cdf((c.study.fovHalfAngleRad/c.study.angleSigmaRad)^2,2, ...
        sum(c.study.angleBiasRad.^2)/c.study.angleSigmaRad^2,'upper');
else
    pOut=double(norm(c.study.angleBiasRad)>c.study.fovHalfAngleRad);
end
m.aber=(1-pOut)*m.aber+.5*pOut; m.op=(1-pOut)*m.op+pOut;
m.meanReceivedPowerW=(1-pOut)*m.meanReceivedPowerW;
m.fovOutageExact=pOut; m.fovOutageSample=1-mean(accepted);
end
