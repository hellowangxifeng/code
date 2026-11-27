function c = rd_config(profile)
%RD_CONFIG All inputs of the dynamic-sea RIS versus rough-sea direct comparison.
%   c = rd_config()          production settings (README section 5)
%   c = rd_config('smoke')   tiny end-to-end check written to results/smoke
%
%   Optics, LC array, device limits, water, beam, phase families, receiver
%   jitter/FOV and the EGG x Gamma-Gamma channel are inherited unchanged from
%   RIS_phase (rp_stage2_config). Only the fields set below are new or
%   changed. Edit the returned structure and pass it to rd_run(c).
if nargin==0, profile = 'production'; end
profile = validatestring(profile,{'production','smoke'});
rd_setup();
c = rp_stage2_config();                 % RIS_phase stage-two model (read-only)
paths = rd_setup();                     % keep RIS_dynamic first on the path
c.version = 'ris-dynamic-v3';
c.profile = profile;
c.paths = paths;
c.output = fullfile(paths.dynamic,'results');
c.cache.enabled = true;
c.cache.folder = fullfile(c.output,'cache');   % never inside RIS_phase

%% Scene: RIS_phase main protocol
c.scene.distancesM = [7.5 10];          % receiver plane above the board top (m)
c.scene.incidenceDeg = 40;              % water-side incidence of the TX axis (deg)

%% Receiver and channel (RIS_phase Fig2/Fig3 FOV; all else from rp_stage2_config)
c.study.fovHalfAngleRad = 0.130;        % FOV half-angle (rad)
c.study.turbulence = 'weak';            % ocean EGG used for design and training
c.study.atmosphere = 'weak';            % atmospheric Gamma-Gamma
c.receiver.minimumCalmPowerW = 1e-6;    % calm ballistic power below this -> pilot alignment
c.receiver.pilotFrames = 12;            % rough-sea frames used for pilot alignment
c.receiver.pilotPhotons = 1e5;

%% Sea states: one directional sea per wind drives the RIS AND the direct surface
c.sea.windsMps = [3 4 5 6];             % wind speeds (m/s), RIS_float range 3-25; calm U = 0 is added
c.sea.windDirectionDeg = 50;            % as W_RIS / RIS_phase conditions
c.ris.windAngleDeg = c.sea.windDirectionDeg;
c.ris.duration = 30;                    % s of six-DOF motion per wind
c.ris.dt = 0.05; c.ris.startTime = 25; c.ris.seed = 600;

%% Dynamic Monte Carlo, identical for both links
c.dynamic.frameStepS = 0.5;             % frozen optical frames (multiple of c.ris.dt)
c.dynamic.photons = 5e5;                % packets per link and frame
c.dynamic.chunkSize = 1e5;
c.dynamic.samples = 2048;               % receiver jitter/angle samples per frame
c.dynamic.calmFrames = 3;               % independent MC repeats of the calm sea
c.grid.bins = [100 100]; c.grid.planeScale = 1;   % wris_mc internal histograms only

%% Rough water-air boundary (direct link)
c.surface.gridStep = 0.0125;            % m, height-field sampling (lambda_min ~ 0.17 m)
c.surface.gridHalfWidth = 3;            % m around the beam crossing; exact Fourier sum outside
c.surface.marchStep = 0.02;             % m, minimum ray-march step
c.surface.rootTolerance = 2e-6;         % m, crossing bisection tolerance
c.surface.crestBlocking = true;         % transmitted rays that re-enter a crest are lost

%% Packet receiver numerics
c.compress.positionM = 2e-3;            % packet merging cell on the receiver plane
c.compress.directionCos = 2e-3;         % packet merging step in direction cosines
c.capture.maxPairs = 2e6;               % sample x packet block size (memory)

%% Sea-aware width selection (trained on an independent sea realization)
c.training.enabled = true;
c.training.objective = 'aber';          % 'aber' or 'op', evaluated at c.study.designPowerW
c.training.families = {'robust','flat_top'};
c.training.seedOffset = 1000;           % wave-phase seed = c.ris.seed + seedOffset
c.training.frameStepS = 1;
c.training.photons = 2e5;
c.training.samples = 1024;

%% Incidence sweep (RIS_phase Fig10 counterpart)
c.incidence.enabled = true;
c.incidence.anglesDeg = [30 35 40 44 47 49 51 54];
c.incidence.windMps = 5;                % must be one of c.sea.windsMps (sea-aware design)
c.incidence.frameStepS = 1;
c.incidence.photons = 3e5;
c.incidence.samples = 1024;

%% Receiver-plane irradiance maps
c.maps.bins = 161;
c.maps.halfWidthDirectM = 2.0;
c.maps.halfWidthRisM = 0.6;
c.maps.snapshotTimeS = 15;

%% Post-processing only (changing these never triggers a recomputation)
c.metrics.turbulence = {'weak','weak';'strong','strong'};  % {ocean, atmosphere}; row 1 = main
c.metrics.fixedPowerDBm = 30;           % summary / wind / incidence figures
c.metrics.cdfEdgesDB = -90:0.5:0;

%% Run control
c.run.resume = true;                    % reuse checkpoints whose inputs and code match
c.run.progressEvery = 5;

if strcmp(profile,'smoke')
    c.scene.distancesM = 7.5;
    c.sea.windsMps = 3;
    c.ris.duration = 2;
    c.mc.photons = 3e4;                 % RIS_phase design-scene MC (rp_mc)
    c.study.widthsM = [.06 .10 .14];
    c.study.centerMargins = [0 .012];
    c.study.samples = 512;
    c.dynamic.frameStepS = 1;
    c.dynamic.photons = 2e4;
    c.dynamic.chunkSize = 1e4;
    c.dynamic.samples = 256;
    c.dynamic.calmFrames = 1;
    c.training.photons = 1e4;
    c.training.samples = 128;
    c.incidence.anglesDeg = [40 50];
    c.incidence.windMps = 3;
    c.incidence.photons = 1e4;
    c.incidence.samples = 128;
    c.receiver.pilotFrames = 2;
    c.receiver.pilotPhotons = 5e3;
    c.surface.gridHalfWidth = 2;
    c.output = fullfile(paths.dynamic,'results','smoke');
    c.cache.folder = fullfile(c.output,'cache');
end
end
