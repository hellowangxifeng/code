function h=rp_aperture_channel(c,g,w,center,s,options)
%RP_APERTURE_CHANNEL Tilted finite detector integral, FOV and Rayleigh jitter.
% I is horizontal flux. Project the physical disk ALONG the nominal ray onto
% the horizontal map; its Jacobian is |n_RX dot s_ray| / s_ray,z.
if nargin<6, options=c.study; end
ray=[center-g.pRis(1:2).',g.rxPlaneZ-g.pRis(3)-c.ris.thickness];
ray=ray/norm(ray);
e1=cross([0 1 0],ray); e1=e1/norm(e1); e2=cross(ray,e1);
jitter=options.positionSigmaM*s.positionNormal;
errors=options.angleBiasRad+options.angleSigmaRad*s.angleNormal;
psi=vecnorm(errors,2,2); tangent=errors./max(psi,realmin);
% The rotation axis s x direction rotates s toward the erroneous normal.
axis=tangent(:,1).*e2-tangent(:,2).*e1;
cosine=cos(psi); sine=sin(psi);
u=cosine.*e1+sine.*cross(axis,repmat(e1,s.n,1),2)+ ...
    (1-cosine).*sum(axis.*e1,2).*axis;
v=cosine.*e2+sine.*cross(axis,repmat(e2,s.n,1),2)+ ...
    (1-cosine).*sum(axis.*e2,2).*axis;
% Disk area quadrature: t=(rho/a)^2 uniform, azimuth uniform.
radialOrder=4; azimuthCount=16;
if isfield(options,'apertureRadialOrder'), radialOrder=options.apertureRadialOrder; end
if isfield(options,'apertureAzimuthCount'), azimuthCount=options.apertureAzimuthCount; end
[nodes,weights]=rp_gauss_legendre(radialOrder);
radius=c.rx.apertureRadius*sqrt((nodes+1)/2); weights=weights/2;
azimuth=2*pi*((1:azimuthCount)-0.5)/azimuthCount;
base=center+jitter; captured=zeros(s.n,1); outside=zeros(s.n,1);
for k=1:numel(radius)
    for j=1:numel(azimuth)
        disk=radius(k)*(cos(azimuth(j))*u+sin(azimuth(j))*v);
        xy=base+disk(:,1:2)-disk(:,3)*ray(1:2)/ray(3);
        intensity=interp2(w.x,w.y,w.irradianceWm2,xy(:,1),xy(:,2),'linear',0);
        captured=captured+intensity*weights(k)/numel(azimuth);
        outside=outside+double(xy(:,1)<w.x(1) | xy(:,1)>w.x(end) | ...
            xy(:,2)<w.y(1) | xy(:,2)>w.y(end))*weights(k)/numel(azimuth);
    end
end
projected=max(0,cosine)/ray(3);
withoutFOV=captured*pi*c.rx.apertureRadius^2.*projected;
accepted=psi<=options.fovHalfAngleRad & cosine>0;
gain=withoutFOV.*accepted/c.source.power;
if options.angleSigmaRad>0
    exactFloor=ncx2cdf((options.fovHalfAngleRad/options.angleSigmaRad)^2,2, ...
        sum(options.angleBiasRad.^2)/options.angleSigmaRad^2,'upper');
else
    exactFloor=double(norm(options.angleBiasRad)>options.fovHalfAngleRad);
end
h=struct('gain',gain,'withoutFOV',withoutFOV/c.source.power, ...
    'fovOutageSample',mean(~accepted),'fovOutageExact',exactFloor, ...
    'nominalDirection',ray,'positionOffsetM',jitter,'angleErrorRad',psi, ...
    'outsideMapApertureFraction',mean(outside),'nominalCenterXY',center);
end
