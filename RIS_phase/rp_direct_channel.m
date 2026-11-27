function h=rp_direct_channel(c,g,w,worldCenter,samples)
%RP_DIRECT_CHANNEL Same physical disk/jitter after horizontal map rotation.
Q=w.worldToMap; rotation=blkdiag(Q,1); mapCenter=worldCenter*Q;
ray=[worldCenter-g.pRis(1:2).',g.rxPlaneZ-g.pRis(3)]; ray=ray/norm(ray);
mapRay=ray*rotation;
e1=cross([0 1 0],ray); e1=e1/norm(e1); e2=cross(ray,e1);
v1=cross([0 1 0],mapRay); v1=v1/norm(v1); v2=cross(mapRay,v1);
angleRotation=([e1;e2]*rotation)*[v1;v2].';
mapped=samples; mapped.positionNormal=samples.positionNormal*Q;
mapped.angleNormal=samples.angleNormal*angleRotation;
opts=c.study; opts.angleBiasRad=c.study.angleBiasRad*angleRotation;
mapGeometry=g; mapGeometry.pRis(1:2)=(g.pRis(1:2).'*Q).';
% The direct interface is the water surface, with no 1 cm glass on top.
mapGeometry.pRis(3)=g.pRis(3)-c.ris.thickness;
h=rp_aperture_channel(c,mapGeometry,w,mapCenter,mapped,opts);
h.nominalDirection=ray; h.nominalCenterXY=worldCenter;
h.positionOffsetM=h.positionOffsetM*Q.';
end
