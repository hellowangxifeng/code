function q = rd_condition(c,s,z)
%RD_CONDITION All five links at one RX height z on one frozen sea (scene s).
%   The receiver plane is Z = Zc + z (Zc: calm LC-screen centre). The RX
%   centre is where the board-local exit direction s.q0 (rd_rx_direction)
%   from the LC centre meets that plane, i.e. inside the cone the RIS can
%   reach. ONE physical receiver - centre, boresight towards the board top
%   centre, and the same Sobol jitter / angle-error samples - is used by
%     1  Direct (no RIS)  HG MC through the frozen rough surface (rd_direct_mc)
%     2  Steering-only    global LC ramp aimed at the RX
%     3+ Flat top / Robust designs over c.design.widthsM, aimed at the RX
%   Links 2-3+ are computed with the complex field and checked with MC
%   (rd_candidate). Scheme selection (Flat top, Robust ABER/OP, Max mean
%   power) is post-processing (rd_evaluate), so it can change freely.
Z = s.Zc+z;
g = s.g;
g.rxPlaneZ = Z;
placement = s.placement;
placement.rxPlaneZ = Z;
[target,ok] = rp_carrier_landing(c,placement,[0 0],s.q0);
assert(ok,'rd:Receiver','The receiver direction misses the plane z = %g m.',z);
g.rxDesired = [target(:);Z];
samples = rp_link_samples(c);
origin = [g.pRis(1:2).',g.pRis(3)+c.ris.thickness];   % = rp_aperture_channel boresight origin
rx = rd_receiver(c,origin,target,Z,samples);
fprintf('[condition] U = %g m/s, z = %g m: RX at (%.3f, %.3f, %.3f) m, %.1f deg off vertical\n', ...
    s.wind,z,target(1),target(2),Z,acosd(rx.ray(3)));

land = rd_to_plane(s.air,Z);
edges = @(x0) x0+linspace(-c.maps.halfWidthM,c.maps.halfWidthM,c.maps.bins+1);
direct = struct('gain',single(rd_capture(c,rd_compress(land,c.compress),rx)), ...
    'map',single(rd_map(land,edges(target(1)),edges(target(2)))), ...
    'mapHalfWidthM',c.maps.halfWidthM,'planePowerW',sum(land.power), ...
    'planeBallisticPowerW',sum(land.power(land.ballistic)),'ledger',s.directLedger);
fprintf('    direct   mean gain %.4g (plane power %.4g W)\n',mean(direct.gain),direct.planePowerW);

specs = struct('family','steering','width',0);
for family = c.design.families
    for width = c.design.widthsM
        specs(end+1) = struct('family',family{1},'width',width); %#ok<AGROW>
    end
end
fieldHash = rd_code_hash('field');
physical = rmfield(rd_physical(c),'direct');   % designs never depend on the no-RIS MC
folder = fullfile(c.cache.folder,'designs');
E = cell(1,numel(specs));
job = @(j) rd_cached(fullfile(folder,sprintf('design_U%g_z%g_%s_w%03d.mat',s.wind,z, ...
    specs(j).family,round(1000*specs(j).width))), ...
    struct('physical',physical,'wind',s.wind,'z',z,'spec',specs(j),'code',fieldHash), ...
    @() rd_candidate(c,s,g,target,specs(j),samples,rx),c.run.resume);
if c.run.parallel && ~isempty(ver('parallel'))
    parfor j = 1:numel(specs)
        E{j} = job(j); %#ok<PFBNS>
    end
else
    for j = 1:numel(specs)
        E{j} = job(j);
    end
end
q = struct('distanceM',z,'wind',s.wind,'Z',Z,'targetXY',target,'rxOrigin',origin, ...
    'rxRay',rx.ray,'psi',single(rx.psi),'fovOutageExact',rx.fovOutageExact, ...
    'pose',s.pose,'q0',s.q0,'direct',direct,'candidates',[E{:}]);
end
