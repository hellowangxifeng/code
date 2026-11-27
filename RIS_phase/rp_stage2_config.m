function c=rp_stage2_config()
%RP_STAGE2_CONFIG Stage-two optical design and joint EGG/Gamma-Gamma channel.
c=rp_config(); c.version='stage2-egg-gg-v1';
c.output=fullfile(c.root,'plots','stage2');
c.cache.folder=fullfile(c.root,'plots','cache','stage2');
c.study.distancesM=[5 7.5 10]; c.study.widthsM=.04:.02:.20;
c.study.families={'gaussian','flat_top','robust'};
c.study.extendedWidthsM=.22:.02:.32;
c.study.jitterExtendedWidthsM=.22:.02:.40;
c.study.lowerWidthsM=[.01 .02 .03];
c.study.minimumWidthM=.04; c.study.maximumAxisRatio=1.35;
c.study.directionSlackDeg=1.5; c.study.centerMargins=[0 .012 .03];
c.study.maximumTailPower=1e-5; c.study.referenceWidthM=.10;
c.study.robustExponent=1.6;
c.study.designPowerW=1; c.study.powerDBm=linspace(-5,45,101);
c.study.positionSigmaM=.080; c.study.angleSigmaRad=.025;
c.study.angleBiasRad=[.025 0]; c.study.fovHalfAngleRad=.240;
c.study.sensitivityFovHalfAngleRad=.090;
c.study.responsivityAW=.5; c.study.noiseStdA=5e-5; c.study.snrThreshold=10;
c.study.turbulence='weak';
% Path-effective Gamma-Gamma shapes, not an invented Cn^2 at a 10 m link.
c.study.atmosphere='weak';
c.study.ggWeak=[17 16]; c.study.ggStrong=[4.2 1.4];
c.study.samples=8192; c.study.validationSamples=32768;
c.study.seed=12112026; c.study.quadratureOrder=32;
c.study.positionSweepM=linspace(.020,.120,26);
c.study.fovSweepRad=linspace(.030,.150,25);
c.study.angleSweepRad=linspace(.010,.060,26);
c.study.jitterGridM=linspace(0,.120,31);
c.study.angleGridDeg=[20 30 40 45 47 48 49 50 52 56 60 65];
c.study.proxyBins=192; c.study.waveBins=128;
c.study.waveCheckBins=192; c.study.raySamplesPerAxis=3;
c.reachability.incidencePolarDeg=40;
az=deg2rad(c.reachability.incidenceAzimuthDeg);
c.geometry.rotationXYDeg=beamRotation(40,az);
c.spot.mapBins=[128 128]; c.mc.photons=300000;
c.wave.maximumOmittedPowerFraction=1e-5;
c.plot.resolution=300; c.plot.tickFontSize=10;
end

function angles=beamRotation(polar,az)
d=[sind(polar)*cos(az);sind(polar)*sin(az);cosd(polar)];
roll=-asind(d(2)); pitch=atan2d(d(1),d(3)); angles=[roll pitch];
end
