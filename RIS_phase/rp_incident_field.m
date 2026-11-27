function f = rp_incident_field(c,g,xy,glassHeight)
%RP_INCIDENT_FIELD PDF Eq.(8) amplitude with an explicitly assumed Gaussian phase.
% |field|^2 is normal power flux in W/m^2, not electric-field V/m squared.
% Gaussian phase is specified by the source, never inferred from an MC histogram.
arguments
    c (1,1) struct
    g (1,1) struct
    xy (:,2) double {mustBeReal,mustBeFinite}
    glassHeight (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = g.hIn
end
nw=c.optics.nWater; ng=c.optics.nGlass; k0=2*pi/c.source.wavelength;
depth=-g.sourceLocal(3);
path=rp_layer_path(xy-g.sourceLocal(1:2).',[depth,glassHeight], ...
    [nw,ng],c.numerics.bisectionIterations);
q=path.q; kg=path.kz(:,2);
bottom=xy-glassHeight*q./kg;
delta=[bottom,zeros(size(bottom,1),1)]-g.sourceLocal.';
axial=delta*g.directionLocal;
radialSquared=max(0,sum(delta.^2,2)-axial.^2);
zR=pi*nw*c.source.waist^2/c.source.wavelength;
widthSquared=c.source.waist^2*(1+(axial/zR).^2);
inverseR=axial./(axial.^2+zR^2);
waterOpl=nw*(axial+radialSquared.*inverseR/2);
phase=k0*(waterOpl+ng*path.length(:,2))-atan2(axial,zR);
irradiance=rp_bsf_ballistic(c,g,xy,glassHeight);
f=struct('field',sqrt(irradiance).*exp(1i*phase),'phaseRad',phase, ...
    'irradianceWm2',irradiance,'bottomXY',bottom,'q',q,'opl',path.opl, ...
    'widthM',sqrt(widthSquared),'rayleighRangeM',zR, ...
    'amplitudeModel','BSF_note Eq.(8), finite-waist 3x3 cubature; no MC fit', ...
    'phaseAssumption','coherent Gaussian wavefront plus glass eikonal');
end
