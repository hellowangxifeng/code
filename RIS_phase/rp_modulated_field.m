function f = rp_modulated_field(c,g,a,p,unitIndex,xy)
%RP_MODULATED_FIELD Expose the ballistic phasor before/after a physical LC tile.
% The glass/air interface and subsequent diffraction are applied by rp_wave.
arguments
    c (1,1) struct
    g (1,1) struct
    a (1,1) struct
    p (1,1) struct
    unitIndex (1,1) double {mustBeInteger,mustBePositive}
    xy (:,2) double {mustBeReal,mustBeFinite}
end
assert(unitIndex<=numel(p.unitIds),'rp:UnitIndex','Unit index is outside the loaded phase array.');
input=rp_incident_field(c,g,xy);
delta=xy-p.referenceXY(unitIndex,:);
k0=2*pi/c.source.wavelength;
profile=rp_phase_profile(c,p,unitIndex,delta);
phase=p.phaseAtReferenceRad(unitIndex)+k0*profile.phaseOverK0;
center=a.unitCenters(p.unitIds(unitIndex),:);
inside=all(abs(xy-center)<=c.array.apertureM/2,2);
incident=input.field*sqrt(a.ballisticModelScale);
after=incident.*exp(1i*phase).*inside;
f=struct('incidentField',incident,'phaseModulationRad',phase,'afterRisField',after, ...
    'irradianceWm2',abs(after).^2,'insideAperture',inside);
end
