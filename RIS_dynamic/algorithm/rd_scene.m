function s = rd_scene(c,wind)
%RD_SCENE One frozen instant of one sea state, shared by every height and link.
%   s = rd_scene(c,wind)
%   1. Calm reference: wris_reference board at its flotation equilibrium,
%      TX 10 m below at c.scene.incidenceDeg (RIS_phase rp_stage2_scene).
%   2. Sea: one directional Neumann sea (sea_wave) drives the six-DOF board
%      motion (RIS_float, rd_motion). The instant t* = c.sea.freezeTimeS
%      is frozen: board pose R(t*), lower-face centre p(t*), and the same
%      wave components give the surface eta(x,y,t*) of the no-RIS link.
%   3. TX: fixed position; with c.tx.trackBoard it re-points its axis at
%      p(t*) (pointing tracking), so the beam stays on the 1 m board.
%   4. RIS optics on the frozen board (RIS_phase, unchanged): lower-face MC
%      (rp_mc), tile illumination (rp_illuminate), tile references
%      (rp_phase). These feed the complex field of every RIS design.
%   5. RX direction q0: centre-ward edge of the cone every tile can reach
%      with its +-15.4 deg commands, plus c.rx.carrierMargin (rd_rx_direction).
%   6. Independent full MC of the RIS input (validation) and the no-RIS
%      HG Monte Carlo through eta(x,y,t*) (rd_direct_mc).
g0 = calmGeometry(c);
m = rd_motion(c,wind,c.ris.seed);
[~,k] = min(abs(m.t-c.sea.freezeTimeS));
pRis = m.bottomCenter(:,k);
R = m.rotationMatrix(:,:,k);
pTx = g0.pTx;
d = g0.dTx;
RTx = g0.RTx;
if c.tx.trackBoard
    d = (pRis-pTx)/norm(pRis-pTx);
    RTx = alignment(g0.dTx,d)*g0.RTx;
end
Zc = g0.pRis(3)+g0.hIn;                 % calm LC-screen centre: heights are measured from here
probeZ = Zc+max(c.scene.distancesM);
motion = struct('t',m.t(k),'rotationMatrix',R,'bottomCenter',pRis);
tx = struct('position',pTx,'direction',d,'rotation',RTx);
rx = struct('position',[pRis(1:2);probeZ],'nominal',[pRis(1:2);probeZ]);
g = rp_geometry(c,motion,tx,rx,1);
fprintf(['[scene] U = %g m/s, t* = %.2f s: tilt %.2f deg, heave %+.3f m, offset (%+.3f, %+.3f) m, ' ...
    'incidence on board %.2f deg\n'],wind,m.t(k),acosd(min(1,R(3,3))),pRis(3)-g0.pRis(3), ...
    pRis(1)-g0.pRis(1),pRis(2)-g0.pRis(2),g.incidenceDeg);

mc = rp_mc(c,g);
a = rp_illuminate(c,g,mc);
base = rp_phase(c,g,a);
if strcmp(c.rx.placement,'frozen')
    placement = g;
    q0 = rd_rx_direction(c,a,base);
else
    mc0 = rp_mc(c,g0);
    a0 = rp_illuminate(c,g0,mc0);
    placement = g0;
    q0 = rd_rx_direction(c,a0,rp_phase(c,g0,a0));
end

fprintf('  RIS validation MC: %d packets\n',round(c.validation.photons));
[~,packets] = wris_mc(rd_transport(c,c.validation.photons),g,c.mc.seed+c.validation.seedOffset);
lcp = rd_ris_packets(c,g,packets);
clear packets

crossing = pTx(1:2)-pTx(3)/d(3)*d(1:2);
sea = rd_surface(m,m.t(k),c.surface,crossing);
fprintf('  no-RIS rough-surface MC: %d packets\n',round(c.direct.photons));
air = rd_direct_mc(rd_transport(c,c.direct.photons),g,sea,c.mc.seed+c.direct.seedOffset);
[eta,slope] = RIS_wave_field(m,crossing.',m.t(k));
if abs(air.energyResidual)>1e-9
    warning('rd:Energy','Direct MC energy ledger residual %.3g at U = %g m/s.',air.energyResidual,wind);
end

pose = struct('timeS',m.t(k),'tiltDeg',acosd(min(1,R(3,3))),'heaveM',pRis(3)-g0.pRis(3), ...
    'offsetXYM',(pRis(1:2)-g0.pRis(1:2)).','normal',R(:,3).','incidenceOnBoardDeg',g.incidenceDeg, ...
    'txRepointDeg',acosd(min(1,dot(d,g0.dTx))),'surfaceEtaM',eta,'surfaceSlopeDeg',atand(norm(slope)), ...
    'significantWaveHeightM',4*sqrt(sum(abs(m.wave.coefficient).^2)/2), ...
    'maxTiltRecordDeg',m.validity.maxTiltDeg,'boardInputPowerW',a.activePowerW, ...
    'boardBallisticPowerW',a.activeDirectPowerW);
ledger = rmfield(air,{'pos','dir','power','ballistic'});
s = struct('wind',wind,'g',g,'g0',g0,'Zc',Zc,'q0',q0,'placement',placement, ...
    'mc',mc,'a',a,'base',base,'lcp',lcp,'air',air,'pose',pose,'directLedger',ledger, ...
    'crossingXY',crossing.');
end

function g = calmGeometry(c)
% rp_study_scene without its Monte Carlo: the motionless reference board.
g = wris_reference(c);
g.sourceLocal = g.RRis.'*(g.pTx-g.pRis);
g.directionLocal = g.RRis.'*g.dTx;
g.hIn = c.ris.thickness*c.optics.screenFraction;
g.hOut = c.ris.thickness-g.hIn;
g.zTop = c.ris.thickness;
g.rxPlaneZ = g.pRis(3)+g.hIn+max(c.scene.distancesM);
g.rxDesired = [g.pRis(1:2);g.rxPlaneZ];
g.time = 0;
g.incidenceDeg = acosd(g.directionLocal(3));
end

function Q = alignment(a,b)
% Smallest rotation taking unit vector a onto unit vector b (Rodrigues).
v = cross(a,b);
s = norm(v);
if s<1e-12
    Q = eye(3);
    return
end
K = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0]/s;
theta = atan2(s,dot(a,b));
Q = eye(3)+sin(theta)*K+(1-cos(theta))*(K*K);
end
