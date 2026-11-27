function [g,mc,a,c]=rp_stage2_scene(c,z,incidence)
%RP_STAGE2_SCENE Explicit incidence at fixed 10 m underwater slant length.
if nargin<3, incidence=c.reachability.incidencePolarDeg; end
az=deg2rad(c.reachability.incidenceAzimuthDeg);
direction=[sind(incidence)*cos(az);sind(incidence)*sin(az);cosd(incidence)];
c.geometry.rotationXYDeg=[-asind(direction(2)),atan2d(direction(1),direction(3))];
c.reachability.incidencePolarDeg=incidence;
[g,mc,a]=rp_study_scene(c,z);
assert(abs(g.incidenceDeg-incidence)<1e-8,'rp:Incidence','Incidence setup mismatch.');
end
