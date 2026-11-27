function r = rp_trace(c,g,a,p)
%RP_TRACE Propagate every illuminated tile sample under its tile's phase ramp.
arguments
    c (1,1) struct
    g (1,1) struct
    a (1,1) struct
    p (1,1) struct
end
m=a.sampleUnitIndex; delta=a.lcXY-p.referenceXY(m,:);
profile=rp_phase_profile(c,p,m,delta);
gradient=profile.gradientOverK0;
q = a.qIn+gradient;
glassSquared = c.optics.nGlass^2-sum(q.^2,2);
glassKz = sqrt(max(glassSquared,realmin));
topXY = a.lcXY+g.hOut*q./glassKz;
inside = abs(topXY(:,1))<=c.ris.Lx/2 & abs(topXY(:,2))<=c.ris.Ly/2;
topWorld = [topXY,repmat(g.zTop,size(q,1),1)]*g.RRis.'+g.pRis.';
airSquared = c.optics.nAir^2-sum(q.^2,2);
direction = [q,sqrt(max(0,airSquared))]/c.optics.nAir*g.RRis.';
up = direction(:,3)>c.numerics.forwardCosTolerance;
distance = (g.rxPlaneZ-topWorld(:,3))./max(direction(:,3),realmin);
valid = glassSquared>0 & airSquared>0 & inside & up & distance>0;
landing = topWorld+distance.*direction; landing(~valid,:) = NaN;
status = ones(size(valid),'uint8'); % 1: no propagating air/glass mode
status(glassSquared>0 & airSquared>0 & inside & (~up | distance<=0)) = 2;
status(~inside) = 3; status(valid) = 0;
r = struct('qOut',q,'loadedGradientOverK0',gradient,'topXY',topXY,'topWorld',topWorld, ...
    'directionWorld',direction,'landing',landing,'valid',valid,'status',status, ...
    'powerW',a.powerW.*rp_fresnel(c.optics.nGlass,c.optics.nAir,q), ...
    'inputPowerW',a.powerW,'unitId',a.unitId);
end
