function result=rp_mc_water_air(c,g,mc)
%RP_MC_WATER_AIR Refract saved lower-face MC packets at a bare sea surface.
arguments
    c (1,1) struct
    g (1,1) struct
    mc (1,1) struct
end
% No glass, tile mask, fill factor, phase or coherent aperture is used here.
q=c.optics.nWater*mc.directionLocal(:,1:2);
T=rp_fresnel(c.optics.nWater,c.optics.nAir,q);
localAir=[q/c.optics.nAir,sqrt(max(0,1-sum(q.^2,2)/c.optics.nAir^2))];
air=localAir*g.RRis.';
surface=[mc.bottomXY,zeros(numel(T),1)]*g.RRis.'+g.pRis.';
travel=(g.rxPlaneZ-surface(:,3))./max(air(:,3),realmin);
valid=T>0 & air(:,3)>1e-12 & travel>0 & isfinite(travel);
xy=surface(valid,1:2)+air(valid,1:2).*reshape(travel(valid),[],1);
power=mc.powerW(valid).*T(valid);
result=struct('landingXY',xy,'powerW',power,'ballistic',mc.ballistic(valid), ...
    'inputPacketIds',find(valid),'directionAir',air(valid,:), ...
    'incidentPowerW',sum(mc.powerW),'reflectedPowerW',sum(mc.powerW.*(1-T)), ...
    'transmittedPowerW',sum(mc.powerW.*T),'receiverPowerW',sum(power), ...
    'nonForwardPowerW',sum(mc.powerW(~valid).*T(~valid)), ...
    'photons',mc.photons,'sourcePowerW',c.source.power, ...
    'model','saved underwater MC packets -> bare water-air Fresnel/Snell -> same world receiver plane');
end
