function m=rp_fig12_metrics(c,h,powerW,order)
%RP_FIG12_METRICS Exact FOV mass with sampled conditional aperture response.
arguments
    c (1,1) struct
    h (1,1) struct
    powerW (1,:) double {mustBePositive}
    order (1,1) double = 64
end
accepted=h.angleErrorRad<=c.study.fovHalfAngleRad & cos(h.angleErrorRad)>0;
floorProbability=h.fovOutageExact;
if ~any(accepted)
    assert(floorProbability>1-1e-12,'rp:ConditionalSamples','No samples inside nonempty FOV.');
    m=struct('aber',.5*ones(size(powerW)),'op',ones(size(powerW)), ...
        'meanReceivedPowerW',zeros(size(powerW)));
    return
end
conditional=struct('gain',h.withoutFOV(accepted),'fovOutageExact',0,'fovOutageSample',0);
m=rp_link_metrics(c,conditional,powerW,order);
m.aber=.5*floorProbability+(1-floorProbability)*m.aber;
m.op=floorProbability+(1-floorProbability)*m.op;
m.meanReceivedPowerW=(1-floorProbability)*m.meanReceivedPowerW;
end
