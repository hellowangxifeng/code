function c = wrd_config(profile)
%WRD_CONFIG Shared inputs for the rough-sea direct link and the floating RIS link.
%   c = wrd_config()             production settings (see README.md)
%   c = wrd_config('smoke')      1 wind, 1 angle, 1 s, 2000 packets: a quick
%                                end-to-end run that writes to results/smoke
%   Edit the returned structure and pass it to run_wris_dynamic(c).
%
%   Everything that defines the sea, the TX, the water and the MC is copied
%   from W_RIS (wris_config) and RIS_phase (rp_config) so that the RIS link
%   is unchanged; only the comparison, receiver and rough-surface fields are new.
if nargin == 0, profile = 'production'; end
profile = validatestring(profile,{'production','smoke'});
folders = wrd_setup();

original = wris_config('production');   % W_RIS physical defaults (read-only)
c = rp_config();                         % RIS optics, LC array, device limits
c.version = 'ris-dynamic-v2';
c.profile = profile;
c.ris = original.ris;                    % 20 s, dt 0.05 s, startTime 25 s, seed 600
c.source = original.source;              % 1 W, 532 nm, 10 mrad, Gaussian waist
c.water = original.water;                % c = 0.151 1/m, HG g = 0.924
c.mc = original.mc;
c.grid = original.grid;
c.grid.planeScale = 1;
c.grid.bins = [40 40];
c.optics.nWater = c.source.mediumRefractiveIndex;
c.sea = original.sea;
c.tx = original.tx;                      % body-fixed TX + OU dead reckoning (W_RIS)

%% Comparison set-up (same sea realization for all links)
c.comparison.conditions = original.comparison.conditions; % [U (m/s), wind direction (deg)]
c.comparison.rotationXYDeg = [0 40; 0 49]; % TX axis: +z about world x, then world y
c.comparison.frameStep = 0.5;            % s, must be a multiple of c.ris.dt
c.comparison.flatBaseline = true;        % also keep the OLD still-water direct link
c.comparison.risWaveBlocking = false;    % true: crests may also block RIS->RX rays
c.mc.photons = 100000;                   % packets per link and frame
c.mc.chunkSize = 25000;

%% One surveyed receiver shared by all links, fixed before either link is run
% x = 8.3 m samples the nominal 40-degree water-air refracted beam at z = 5 m.
c.receiver.center = [8.3 0 5];
c.receiver.normal = c.receiver.center/norm(c.receiver.center); % boresight = arrival direction
c.receiver.radius = 0.05;                % m, aperture radius
c.receiver.fovDeg = 15;                  % FOV half-angle
c.receiver.mapEdges = linspace(-6,6,241);% detector-plane irradiance map (m)

%% Rough water-air boundary numerics (direct link)
c.surface.gridStep = 0.075;              % m, height-field sampling for interpolation
c.surface.gridHalfWidth = 40;            % m, exact Fourier sum is used outside
c.surface.marchStep = 0.04;              % m, minimum ray-march step
c.surface.rootTolerance = 2e-6;          % m, bisection tolerance of the crossing

%% Statistics over the frames
c.metrics.outageThresholdW = 1e-6;       % W (-30 dBm); set to the real RX sensitivity
c.metrics.displayFloorW = 1e-9;          % W, where zero-hit frames are drawn

%% Run control
c.output = fullfile(folders.dynamic,'results');
c.run.resume = true;                     % reuse finished wind conditions after a crash
c.run.checkpointFolder = '';             % '' -> <output>/checkpoints
c.run.progressEvery = 5;                 % print every n frames

if strcmp(profile,'smoke')
    c.comparison.conditions = c.comparison.conditions(1,:);
    c.comparison.rotationXYDeg = c.comparison.rotationXYDeg(1,:);
    c.ris.duration = 1;
    c.mc.photons = 2000;
    c.mc.chunkSize = 1000;
    c.output = fullfile(folders.dynamic,'results','smoke');
end
end
