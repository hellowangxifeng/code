function c = rp_config()
%RP_CONFIG Physical, discretisation and plot parameters for finite RIS beams.
rp_setup();
base = wris_config('production');
c.version = 'world-z-electrode-1um-v7';
c.root = fileparts(mfilename('fullpath'));
c.ris = base.ris; c.tx = base.tx; c.sea = base.sea;
c.geometry.distance = 10;
c.geometry.rotationXYDeg = [45 45];
c.ris.Lx = 1; c.ris.Ly = 1; c.ris.thickness = 0.01;
c.ris.duration = 20; c.ris.dt = 0.05; c.ris.startTime = 25;
c.ris.seed = 600; c.ris.plot = false;
c.conditions = [3 50;5 50];
c.snapshotTimes = [0 10 20];
c.tx.seed = 111; c.tx.initialYaw = pi/2;
c.tx.meanUVR = [0.008;0;deg2rad(0.03)];
c.tx.stdUVR = [0.012;0.008;deg2rad(0.12)];
c.tx.correlationTime = [2;2;2]; c.tx.pointing = 'body-fixed';
c.optics.nWater = 1.33; c.optics.nGlass = 1.50; c.optics.nAir = 1;
c.optics.screenFraction = 0.5;
c.source.wavelength = 532e-9;
c.source.power = 1;                 % ideal direct-beam power, W
c.source.divergence = 0.01;         % Gaussian 1/e^2 slope radius = tan(divergence)
c.source.waist = c.source.wavelength/(pi*c.optics.nWater*c.source.divergence);
c.water = base.water;
c.mc = base.mc;
c.mc.photons = 300000;
c.mc.chunkSize = 100000;
c.mc.mapBins = [180 180];
c.source.edgeRelativeIntensity = exp(-2); % display contour only; NOT a physical stop
c.array.count = [40 40];            % [x y] tile counts
c.array.pitchM = [0.025 0.025];
c.array.elementSizeM = [0.015 0.015];
c.array.apertureM = [0.0144 0.0144]; % actual transparent LC area per tile
c.device = wris_device(c.source.wavelength); % 1 um = 0.7 um electrode + 0.3 um gap
c.phase.model = 'linear';          % closest realizable world +z carrier; no beam expansion
c.phase.separable = false;         % diagnostic curvature cross term; linear ignores this
c.input.apertureOrder = 6;         % Gauss-Legendre tile-power quadrature
c.reachability.incidencePolarDeg = 60; % fixed WATER incidence, relative to local +z
c.reachability.incidenceAzimuthDeg = atan2d(-sqrt(2),1);
c.reachability.sStep = 0.0025;       % uniform step in air direction cosines sX,sY
c.reachability.detailStep = 3e-4;    % finer step inside the detail window
c.reachability.detailPadding = 0.30; % margin around the command-feasible directions
c.reachability.maxPointingErrorDeg = 0.05; % 3-D angular separation from target
c.rx.height = 5;                    % WORLD horizontal plane above initial top
c.rx.sigma = 0.04; c.rx.correlationTime = 1;
c.rx.seed = 6101; c.rx.startAtNominal = true;
c.rx.apertureRadius = 0.05;
c.spot.minIrradianceWm2 = 0.005;    % measurement only; never changes phase commands
c.spot.coverageFraction = 0.90;
c.spot.mapBins = [220 220];
c.wave.apertureSamples = 96;
c.wave.apertureOversampling = 1;
c.wave.fftPadding = 2;
c.wave.receiverSamplingMargin = 2; % Fourier deposition step <= half a detector pixel
c.wave.maximumFrequencySamples = 32000000; % local CZT spectrum, no dense padding cap
c.wave.maxApertureSamples = 8192;    % per eigen-axis; anisotropic FFT avoids square oversampling
c.wave.supportRelativeFloor = 1e-8;
c.wave.pixelCoherence = true;       % analytic sinc averaging of pair interference
c.wave.minimumHalfWidthM = 0.002;
c.wave.windowPowerFraction = 0.995; % deterministic transmitted ballistic ray envelope
c.wave.windowPadding = 0.15;
c.wave.windowSamples = 5;          % midpoint samples per tile axis, independent of MC hits
c.wave.maximumOmittedPowerFraction = 1e-7; % explicitly account for negligible coherent tails
c.detector.minimumRadiusPixels = 4;
c.detector.maxRefinements = 2;
c.detector.zoomFactor = 0.25;
c.detector.maximumUnresolvedPowerFraction = 1e-3;
c.numerics.bisectionIterations = 48;
c.numerics.forwardCosTolerance = 1e-10;
c.cache.enabled = true; c.cache.folder = fullfile(c.root,'cache');
c.output = fullfile(c.root,'results','validated_optics');
c.plot.visible = 'off'; c.plot.resolution = 180;
c.plot.figureSize = [40 40 1680 900];
c.plot.tickFontSize = 12; c.plot.labelFontSize = 14;
c.plot.annotationFontSize = 9;
c.plot.titleFontSize = 12; c.plot.dynamicRange = 1e4;
end
