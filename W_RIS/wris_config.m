function c = wris_config(profile)
%WRIS_CONFIG Physical defaults; edit run_wris_irradiance.m for a new study.
if nargin==0, profile='production'; end
c.version='wris-8-electrode-1um';
c.root=fileparts(mfilename('fullpath'));
c.geometry.distance=10;
c.geometry.rotationXYDeg=[0 60]; % +z, active world-x then world-y rotation
c.ris=get_default_RIS_params();
c.ris.duration=20; c.ris.dt=0.05; c.ris.startTime=25;
c.ris.fMinHz=0.02; c.ris.fMaxHz=3;
c.ris.frequencyBins=200; c.ris.directionBins=36;
c.ris.seed=600; c.ris.plot=false;
c.comparison.rotationXYDeg=[0 30;0 60;0 80;45 45;-60 30;30 80];
c.comparison.conditions=[3 50;5 50]; % [wind speed (m/s), direction (deg)]
c.comparison.site='center'; % center, crest, trough, rising, falling
c.sampling.comparisonTime=20; % shared clock; appreciable calculated yaw, not amplified
c.sea.quadrature='trapezoid'; c.sea.siteSamples=801;
seaDefaults=get_default_sea_params();
c.sea.neumannConstant=seaDefaults.neumannConstant;
c.tx.seed=111; c.tx.initialYaw=pi/2;
c.tx.meanUVR=[0.008;0;deg2rad(0.03)];
c.tx.stdUVR=[0.012;0.008;deg2rad(0.12)];
c.tx.correlationTime=[2;2;2];
c.tx.pointing='body-fixed';
c.source.power=1; c.source.wavelength=532e-9;
c.device=wris_device(c.source.wavelength);
c.source.mediumRefractiveIndex=1.33;
c.source.divergence=0.01;
c.source.waist=c.source.wavelength/(pi*c.source.mediumRefractiveIndex*c.source.divergence);
c.water.type='c0.151';
[c.water.absorption,c.water.scattering,c.water.attenuation]=Table_abc(c.water.type);
c.water.anisotropy=0.924;
c.turbulence.kind='weak'; c.turbulence.seed=260617;
c.mc.photons=1e6; c.mc.chunkSize=100000;
c.mc.maxScatter=1000; c.mc.maxPath=200;
c.mc.rouletteThreshold=1e-4; c.mc.rouletteSurvival=0.1;
c.mc.seed=20261014;
c.grid.bins=round([c.ris.Lx c.ris.Ly]/0.01); c.grid.planeScale=6;
c.validation.photons=[10000 100000 1000000];
c.validation.repeats=4; c.validation.seed=9500;
c.output=fullfile(c.root,'results','incidence_xy');
if strcmp(profile,'smoke'), c.output=[c.output '_smoke']; end
c.cache.folder=fullfile(c.root,'cache'); c.cache.enabled=true;
c.plot.resolution=300; c.plot.dynamicRange=1e4; c.plot.padding=0.15;
c.plot.rowsPerFigure=3; c.plot.tickFontSize=16; c.plot.labelFontSize=17;
c.plot.titleFontSize=12; c.plot.figureTitleFontSize=16;
c.plot.visible='off'; c.plot.quantity='noTurbulence';
if strcmp(profile,'smoke')
    c.mc.photons=10000;
    c.validation.photons=[2000 10000 30000]; c.validation.repeats=3;
elseif ~strcmp(profile,'production')
    error('wris:Profile','Use production or smoke.');
end
end
