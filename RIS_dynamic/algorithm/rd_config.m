function c = rd_config(profile)
%RD_CONFIG All inputs of the frozen-sea RIS versus no-RIS comparison.
%   c = rd_config()          production settings (README / 算法与图件说明.md)
%   c = rd_config('smoke')   tiny end-to-end check written to results/smoke
%
%   Optics, LC array, device limits, water, beam, phase families, receiver
%   jitter / angle error and the EGG x Gamma-Gamma channel are inherited
%   unchanged from RIS_phase (rp_stage2_config). Only the fields set below
%   are new or changed. Edit the returned structure and pass it to rd_run(c).
if nargin==0, profile = 'production'; end
profile = validatestring(profile,{'production','smoke'});
rd_setup();
c = rp_stage2_config();                 % RIS_phase stage-two model (read-only)
paths = rd_setup();                     % keep RIS_dynamic first on the path
c.version = 'ris-dynamic-v4';
c.profile = profile;
c.paths = paths;
c.output = fullfile(paths.dynamic,'results');
c.cache.enabled = true;
c.cache.folder = fullfile(c.output,'cache');   % never inside RIS_phase

%% Geometry: TX 10 m under the calm board, RX plane z above the calm RIS centre
c.scene.incidenceDeg = 40;              % calm-sea water-side incidence of the TX axis (deg)
c.scene.distancesM = [10 20 30 50];     % RX plane height above the calm LC-screen centre (m)
c.main.distanceM = 30;                  % Fig1-3, Fig6: one height ...
c.main.windMps = 9;                     % ... and one sea state of the grid below
az = deg2rad(c.reachability.incidenceAzimuthDeg);
beam = [sind(c.scene.incidenceDeg)*cos(az);sind(c.scene.incidenceDeg)*sin(az); ...
    cosd(c.scene.incidenceDeg)];
c.geometry.rotationXYDeg = [-asind(beam(2)),atan2d(beam(1),beam(3))];   % = rp_stage2_scene
c.reachability.incidencePolarDeg = c.scene.incidenceDeg;
c.tx.trackBoard = true;                 % TX stays put but re-points at the moving lower face

%% Sea: one directional sea per wind; ONE frozen instant drives board and surface
c.sea.windsMps = [3 9 12];              % RIS_float accepts 3-25 m/s
c.sea.windDirectionDeg = 50;
c.sea.freezeTimeS = 15;                 % frozen optical instant inside the motion record (s)
c.ris.windAngleDeg = c.sea.windDirectionDeg;
c.ris.duration = 20; c.ris.dt = 0.05; c.ris.startTime = 25; c.ris.seed = 600;

%% Receiver: one physical RX shared by all five links
c.rx.placement = 'frozen';              % 'frozen': inside the RIS-reachable cone of the frozen
                                        % board; 'calm': cone of the motionless board
c.rx.carrierMargin = 0.03;              % move the RX into the reachable cone (direction cosines)
c.study.fovHalfAngleRad = 0.130;        % RIS_phase Fig2/3 FOV
c.study.samples = 8192;                 % jitter / angle-error Sobol samples (paired, all links)
c.study.turbulence = 'weak';            % overwritten per row of c.metrics.turbulence
c.study.atmosphere = 'weak';

%% RIS phase designs (RIS_phase families, aimed at the common RX)
c.design.families = {'flat_top','robust'};
c.design.widthsM = [.04 .06 .08 .10 .13 .16 .20 .25 .30];   % target spot widths (m)
% Spot feasibility (RIS_phase): c.study.minimumWidthM, c.study.maximumAxisRatio, r50 >= 2 cm.

%% Complex field: rp_wave per tile, beamlet spectral windows, neighbour coherence
c.wave.receiverHalfWidthM = 0.55;       % map covers RX +- 5 sigma_p jitter + aperture
c.wave.pixelsPerRadius = 4;             % pixel <= aperture radius / 4
c.wave.minBins = 128; c.wave.maxBins = 256;
c.wave.beamletGuardLobes = 64;          % spectral guard around each tile beamlet (lobes of lambda/d)
c.wave.coherentNeighbours = 2;          % exact pixel-averaged cross terms up to this tile ring
c.mc.photons = 3e5;                     % rp_mc packets: scattered halo of the complex field

%% Monte Carlo of the RIS links (independent check of the complex field)
c.validation.photons = 1e6;
c.validation.seedOffset = 7919;
c.grid.bins = [100 100]; c.grid.planeScale = 1;   % wris_mc internal histograms only

%% No-RIS direct link: HG Monte Carlo through the frozen rough surface only
c.direct.photons = 2e6;
c.direct.seedOffset = 31337;
c.dynamic.chunkSize = 1e5;              % packets per MC chunk (all Monte Carlo runs)
c.surface.gridStep = 0.0125;            % m, height-field sampling
c.surface.gridHalfWidth = 3;            % m around the beam crossing; exact Fourier sum outside
c.surface.marchStep = 0.02;             % m, minimum ray-march step
c.surface.rootTolerance = 2e-6;         % m, crossing bisection tolerance
c.surface.crestBlocking = true;         % transmitted rays that re-enter a crest are lost
c.compress.positionM = 2e-3;            % packet merging cell on the receiver plane
c.compress.directionCos = 2e-3;         % packet merging step in direction cosines
c.capture.maxPairs = 2e6;               % sample x packet block size (memory)
c.maps.halfWidthM = 0.6;                % direct-link map around the RX
c.maps.bins = 161;

%% Post-processing only (changing these never triggers a recomputation)
c.metrics.powerDBm = -10:0.5:90;
c.metrics.turbulence = {'weak','weak';'strong','weak';'weak','strong';'strong','strong'};
                                        % {ocean EGG, atmosphere GG}; row 1 = main channel
c.metrics.order = 64;                   % rp_fading_response log-grid refinement
c.metrics.policy = 'per_power';         % RIS_phase Fig2/3: reselect the phase at each power;
                                        % 'design_power': freeze the choice at designPowerDBm
c.metrics.designPowerDBm = 30;          % 1 W: maps, validation table, frozen policy
c.rank.targets = [1e-3 1e-3];           % ABER, OP targets for the required-power ranking
c.rank.toleranceDB = 0.5;               % optimizers within this of the best are co-winners
c.rank.marginDB = 1;                    % "clearly better": ahead of every other link by this

%% Run control
c.run.resume = true;                    % reuse checkpoints whose inputs and code match
c.run.parallel = false;                 % parfor over the designs of one condition

if strcmp(profile,'smoke')
    c.scene.distancesM = [10 30];
    c.sea.windsMps = 3;
    c.main.windMps = 3;
    c.ris.duration = 16;
    c.design.widthsM = [.06 .12 .20];
    c.mc.photons = 3e4;
    c.validation.photons = 5e4;
    c.direct.photons = 5e4;
    c.dynamic.chunkSize = 2.5e4;
    c.study.samples = 512;
    c.wave.maxBins = 128;
    c.surface.gridHalfWidth = 2;
    c.metrics.order = 32;
    c.output = fullfile(paths.dynamic,'results','smoke');
    c.cache.folder = fullfile(c.output,'cache');
end
end
