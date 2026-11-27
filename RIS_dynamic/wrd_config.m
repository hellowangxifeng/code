function c = wrd_config(profile)
%WRD_CONFIG W_RIS inputs and the common physical receiver for both links.
if nargin == 0, profile = 'production'; end
profile = validatestring(profile, {'production','smoke'});
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root,'RIS_phase'),fullfile(root,'W_RIS'), ...
    fullfile(root,'sea_wave'),fullfile(root,'yimiao','fig3-7'));

original = wris_config();
c = rp_config();
c.ris = original.ris;
c.source = original.source;
c.water = original.water;
c.mc = original.mc;
c.grid = original.grid;
c.grid.planeScale = 1;
c.grid.bins = [40 40];
c.optics.nWater = c.source.mediumRefractiveIndex;
c.sea = original.sea;
c.tx = original.tx;
c.comparison.conditions = original.comparison.conditions;
c.comparison.rotationXYDeg = [0 40; 0 49];
c.comparison.frameStep = 0.5;
c.mc.photons = 100000;
c.mc.chunkSize = 25000;

% One surveyed location and boresight, fixed before evaluating either link.
% x=8.3 m samples the nominal 40-degree water-air refracted beam at z=5 m.
c.receiver.center = [8.3 0 5];
c.receiver.normal = c.receiver.center/norm(c.receiver.center); % faces arrival direction
c.receiver.radius = 0.05;
c.receiver.fovDeg = 15; % acceptance half-angle
c.receiver.mapEdges = linspace(-6,6,241);
c.surface.gridStep = 0.075;
c.surface.gridHalfWidth = 40;
c.surface.marchStep = 0.04;
c.surface.rootTolerance = 2e-6;
c.output = fullfile(here,'results');
if strcmp(profile,'smoke')
    c.comparison.conditions = c.comparison.conditions(1,:);
    c.comparison.rotationXYDeg = c.comparison.rotationXYDeg(1,:);
    c.ris.duration = 1;
    c.mc.photons = 2000;
    c.mc.chunkSize = 1000;
    c.output = fullfile(here,'smoke');
end
end
