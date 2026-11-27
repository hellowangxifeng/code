function B = rd_condition(c,wind,S,opts)
%RD_CONDITION Frame-by-frame paired Monte Carlo of all links for one sea state.
%   B = rd_condition(c,wind,S,opts), S from RD_LINKS, opts fields:
%     frameStepS, photons, seedOffset, calmFrames, maps, snapshotTimeS, label.
%   wind = 0 is the calm sea: motionless board, flat water, opts.calmFrames
%   independent MC repeats. For wind > 0, at every frame t_f:
%     1. the board pose comes from the six-DOF motion (RD_MOTION);
%     2. the SAME wave components give the frozen rough surface eta(x,y,t_f);
%     3. RIS link : wris_mc to the moving lower face -> LC masks -> frozen
%                   tile phases (rp_trace) -> receiver plane;
%        direct   : rd_direct_mc through eta(x,y,t_f) -> receiver plane;
%        both links share TX, beam, water and the MC seed of the frame;
%     4. every link is scored by its fixed receiver (rd_capture), giving
%        samples x frames instantaneous gains for RD_METRICS.
transport = rd_transport(c,opts.photons);
if wind==0
    m = rd_static_motion(S.g0,opts.calmFrames);
    keys = 1:numel(m.t);
    Hs = 0;
    validity = struct();
else
    m = rd_motion(c,wind,c.ris.seed);
    keys = rd_frame_keys(c,m,opts.frameStepS);
    Hs = 4*sqrt(sum(abs(m.wave.coefficient).^2)/2);
    validity = m.validity;
end
[tx,rxPose] = rd_poses(S.g0,numel(m.t));
flat = rd_flat_surface();
[nz,K] = size(S.links);
nf = numel(keys);
ns = size(S.links{1,1}.rx.offset,1);
gain = nan(ns,nf,nz,K,'single');
planePowerW = nan(nf,nz,K);
centroidXY = nan(nf,2,nz,K);
bins = c.maps.bins;
if opts.maps
    meanMap = zeros(bins,bins,nz,K,'single');
else
    meanMap = zeros(0,0,'single');
end
snapshotMap = meanMap;
halfWidth = zeros(nz,K);
receivers = cell(nz,K);
for iz = 1:nz
    for k = 1:K
        L = S.links{iz,k};
        halfWidth(iz,k) = L.mapHalfWidthM;
        receivers{iz,k} = struct('key',L.key,'available',L.available,'center',L.rx.center, ...
            'ray',L.rx.ray,'valid',L.rx.valid,'fovOutageExact',L.rx.fovOutageExact, ...
            'fovOutageSample',L.rx.fovOutageSample);
    end
end
snapshot = 1;
if wind>0, [~,snapshot] = min(abs(m.t(keys)-opts.snapshotTimeS)); end
F = struct('t',reshape(m.t(keys),[],1),'tiltDeg',zeros(nf,1),'heaveM',zeros(nf,1), ...
    'slopeDeg',zeros(nf,1),'tirFraction',zeros(nf,1),'ballisticTirFraction',zeros(nf,1), ...
    'ballisticTransmittance',zeros(nf,1),'meanBallisticTiltDeg',zeros(nf,1), ...
    'directTransmittedW',zeros(nf,1),'crestBlockedW',zeros(nf,1), ...
    'risPlatePowerW',zeros(nf,1),'risLcPowerW',zeros(nf,1),'energyResidual',zeros(nf,2));
clock = tic;
for f = 1:nf
    k = keys(f);
    g = rp_geometry(c,m,tx,rxPose,k);
    if wind==0
        sea = flat;
    else
        sea = rd_surface(m,m.t(k),c.surface,S.crossXY);
        [~,slope] = RIS_wave_field(m,S.crossXY,m.t(k));
        F.slopeDeg(f) = atand(norm(slope));
    end
    seed = c.mc.seed+opts.seedOffset+f;
    [plate,packets] = wris_mc(transport,g,seed);
    lcp = rd_ris_packets(c,g,packets);
    air = rd_direct_mc(transport,g,sea,seed);
    F.tiltDeg(f) = acosd(min(1,max(-1,g.RRis(3,3))));
    F.heaveM(f) = g.pRis(3)-S.g0.pRis(3);
    F.tirFraction(f) = air.tirFraction;
    F.ballisticTirFraction(f) = air.ballisticTirFraction;
    F.ballisticTransmittance(f) = air.ballisticTransmittance;
    F.meanBallisticTiltDeg(f) = air.meanBallisticTiltDeg;
    F.directTransmittedW(f) = air.transmittedPowerW;
    F.crestBlockedW(f) = air.crestBlockedPowerW;
    F.risPlatePowerW(f) = plate.capturedPower;
    F.risLcPowerW(f) = sum(lcp.powerW);
    F.energyResidual(f,:) = [plate.energyResidual,air.energyResidual];
    for iz = 1:nz
        direct = rd_to_plane(air,S.Z(iz));
        for kk = 1:K
            L = S.links{iz,kk};
            if ~L.available, continue; end
            if isempty(L.p)
                land = direct;
            else
                land = rd_trace(c,g,lcp,L.p,S.Z(iz));
            end
            planePowerW(f,iz,kk) = sum(land.power);
            spot = rd_spot(land);
            centroidXY(f,:,iz,kk) = spot.centerXY;
            gain(:,f,iz,kk) = single(rd_capture(c,rd_compress(land,c.compress),L.rx));
            if opts.maps && L.rx.valid
                I = single(rd_map(land,L.rx.center,L.mapHalfWidthM,bins));
                meanMap(:,:,iz,kk) = meanMap(:,:,iz,kk)+I/nf;
                if f==snapshot, snapshotMap(:,:,iz,kk) = I; end
            end
        end
    end
    if f==1 || f==nf || mod(f,c.run.progressEvery)==0
        fprintf('    [%s] U = %g m/s | frame %d/%d (t = %5.2f s) | %.1f min\n', ...
            opts.label,wind,f,nf,m.t(k),toc(clock)/60);
    end
end
residual = max(abs(F.energyResidual(:)));
if residual>1e-9
    warning('rd:Energy','MC energy ledger residual %.3g at U = %g m/s.',residual,wind);
end
B = struct('wind',wind,'Hs',Hs,'validity',validity,'frames',F,'gain',gain, ...
    'planePowerW',planePowerW,'centroidXY',centroidXY,'meanMap',meanMap, ...
    'snapshotMap',snapshotMap,'snapshotTimeS',m.t(keys(snapshot)), ...
    'mapHalfWidthM',halfWidth,'receivers',{receivers},'specs',{S.signature}, ...
    'Z',S.Z,'distancesM',S.z,'label',opts.label,'elapsedSeconds',toc(clock));
end
