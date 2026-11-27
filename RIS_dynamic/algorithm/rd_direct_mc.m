function air = rd_direct_mc(c,g,sea,seed)
%RD_DIRECT_MC No-RIS link: HG Monte Carlo through a frozen (rough) sea surface.
%   air = rd_direct_mc(c,g,sea,seed) launches c.mc.photons packets from the
%   TX pose in g (pTx, RTx, dTx) with the source law of W_RIS/wris_mc and
%   transports them through the water column (implicit capture, HG phase
%   function, Russian roulette, source-plane back escape, path/scatter caps)
%   until their FIRST upward crossing of the boundary z = eta(x,y,t) given by
%   sea (RD_SURFACE: the frozen sea of RD_SCENE; a zero-amplitude sea is still water).
%
%   At the crossing point the local normal n = [-grad eta, 1]/|.| gives
%     q = n_w |d - (d.n) n|          (tangential wave vector / k0)
%     d_air = (n_w (d-(d.n)n) + sqrt(n_a^2-q^2) n) / n_a     (vector Snell)
%     T = unpolarised Fresnel transmittance (rp_fresnel), T = 0 for q >= n_a.
%   Reflected and totally reflected power is terminated and logged. With
%   c.surface.crestBlocking, transmitted rays that re-enter the water before
%   leaving the wave slab |z| <= max|eta| are removed (logged as blocked).
%
%   Output air: exit points pos, unit air directions dir, power (W) and the
%   ballistic flag of every transmitted packet, plus the energy ledger.
old = rng;
restore = onCleanup(@() rng(old));
rng(seed,'twister');
photons = c.mc.photons;
absorption = c.water.absorption;
scattering = c.water.scattering;
extinction = absorption+scattering;
L = struct('incident',0,'fresnelReflected',0,'tir',0,'transmitted',0,'absorbed',0, ...
    'backEscape',0,'truncated',0,'rouletteDelta',0,'crestBlocked',0, ...
    'incidentBallistic',0,'tirBallistic',0,'transmittedBallistic',0,'tiltSumBallistic',0);
nChunk = ceil(photons/c.mc.chunkSize);
store = repmat({zeros(0,8)},nChunk,1);
for first = 1:c.mc.chunkSize:photons
    chunk = ceil(first/c.mc.chunkSize);
    n = min(c.mc.chunkSize,photons-first+1);
    uv = c.source.waist/2*randn(n,2);
    slope = c.source.divergence/2*randn(n,2);
    launch = [slope,ones(n,1)];
    direction = (launch./vecnorm(launch,2,2))*g.RTx.';
    position = [uv,zeros(n,1)]*g.RTx.'+g.pTx.';
    weight = ones(n,1);
    path = zeros(n,1);
    collisions = zeros(n,1);
    active = true(n,1);
    outPos = zeros(n,3);
    outDir = zeros(n,3);
    outW = zeros(n,1);
    near = position(:,3)>=-sea.amplitudeBound;
    if any(near)
        assert(all(position(near,3)<sea.evaluate(position(near,1:2))), ...
            'rd:SourceMedium','The TX launch plane must start under water.');
    end
    while any(active)
        ids = find(active);
        p = position(ids,:);
        d = direction(ids,:);
        free = inf(numel(ids),1);
        if extinction>0
            free = -log(max(rand(numel(ids),1),realmin))/extinction;
        end
        remaining = c.mc.maxPath-path(ids);
        axial = d*g.dTx;
        back = ((g.pTx.'-p)*g.dTx)./axial;
        back(axial>=-1e-12 | back<1e-10) = Inf;
        limit = min([free,remaining,back],[],2);
        s = rd_crossing(p,d,limit,sea,c.surface,1);
        hit = isfinite(s);
        reached = ids(hit);
        if any(hit)
            exitPoint = p(hit,:)+d(hit,:).*s(hit);
            [~,gradient] = sea.evaluate(exitPoint(:,1:2));
            normal = [-gradient,ones(nnz(hit),1)];
            normal = normal./vecnorm(normal,2,2);
            [airDirection,T,tir] = refract(c,d(hit,:),normal);
            w = weight(reached);
            ballistic = collisions(reached)==0;
            outPos(reached,:) = exitPoint;
            outDir(reached,:) = airDirection;
            outW(reached) = w.*T;
            L.incident = L.incident+sum(w);
            L.tir = L.tir+sum(w(tir));
            L.transmitted = L.transmitted+sum(w.*T);
            L.fresnelReflected = L.fresnelReflected+sum(w(~tir).*(1-T(~tir)));
            L.incidentBallistic = L.incidentBallistic+sum(w(ballistic));
            L.tirBallistic = L.tirBallistic+sum(w(ballistic & tir));
            L.transmittedBallistic = L.transmittedBallistic+sum(w(ballistic).*T(ballistic));
            L.tiltSumBallistic = L.tiltSumBallistic+ ...
                sum(w(ballistic).*acosd(min(1,normal(ballistic,3))));
        end
        escaped = ~hit & back<=free & back<=remaining;
        capped = ~hit & ~escaped & remaining<=free;
        L.backEscape = L.backEscape+sum(weight(ids(escaped)));
        L.truncated = L.truncated+sum(weight(ids(capped)));
        active(ids(hit | escaped | capped)) = false;
        move = ~hit & ~escaped & ~capped;
        collide = ids(move);
        if isempty(collide), continue; end
        position(collide,:) = p(move,:)+d(move,:).*free(move);
        path(collide) = path(collide)+free(move);
        collisions(collide) = collisions(collide)+1;
        L.absorbed = L.absorbed+sum(weight(collide))*absorption/extinction;
        weight(collide) = weight(collide)*scattering/extinction;
        low = collide(weight(collide)>0 & weight(collide)<c.mc.rouletteThreshold);
        before = weight(low);
        survives = rand(numel(low),1)<c.mc.rouletteSurvival;
        weight(low) = before.*survives/c.mc.rouletteSurvival;
        L.rouletteDelta = L.rouletteDelta+sum(weight(low)-before);
        active(collide(weight(collide)==0)) = false;
        killed = collide(active(collide) & collisions(collide)>=c.mc.maxScatter);
        L.truncated = L.truncated+sum(weight(killed));
        active(killed) = false;
        scatter = collide(active(collide));
        direction(scatter,:) = scatterHG(direction(scatter,:),c.water.anisotropy);
    end
    sent = outW>0;
    P = outPos(sent,:);
    D = outDir(sent,:);
    W = outW(sent);
    B = collisions(sent)==0;
    if c.surface.crestBlocking && sea.amplitudeBound>0 && any(sent)
        blocked = crestBlocked(c,P,D,sea);
        L.crestBlocked = L.crestBlocked+sum(W(blocked));
        P = P(~blocked,:); D = D(~blocked,:); W = W(~blocked); B = B(~blocked);
    end
    store{chunk} = [P,D,W,double(B)];
end
A = vertcat(store{:});
scale = c.source.power/photons;
air = struct('pos',A(:,1:3),'dir',A(:,4:6),'power',A(:,7)*scale,'ballistic',A(:,8)>0, ...
    'incidentPowerW',L.incident*scale,'transmittedPowerW',L.transmitted*scale, ...
    'tirPowerW',L.tir*scale,'crestBlockedPowerW',L.crestBlocked*scale, ...
    'tirFraction',L.tir/max(L.incident,realmin), ...
    'ballisticTirFraction',L.tirBallistic/max(L.incidentBallistic,realmin), ...
    'ballisticTransmittance',L.transmittedBallistic/max(L.incidentBallistic,realmin), ...
    'meanBallisticTiltDeg',L.tiltSumBallistic/max(L.incidentBallistic,realmin), ...
    'energyResidual',(L.transmitted+L.fresnelReflected+L.tir+L.absorbed+L.backEscape+ ...
    L.truncated-photons-L.rouletteDelta)/photons, ...
    'ledger',L,'photons',photons);
end

function blocked = crestBlocked(c,p,d,sea)
% A transmitted ray is lost if it meets the water again before it climbs out
% of the wave slab. Rays that cannot climb (d_z <= 0) are always lost.
blocked = d(:,3)<=1e-9;
ids = find(~blocked);
if isempty(ids), return; end
offset = 8*c.surface.rootTolerance;
start = p(ids,:)+offset*d(ids,:);
limit = min(50,max(0,(sea.amplitudeBound+1e-3-start(:,3))./d(ids,3)));
again = rd_crossing(start,d(ids,:),limit,sea,c.surface,-1);
blocked(ids(isfinite(again))) = true;
end

function [direction,T,tir] = refract(c,incident,normal)
incidence = sum(incident.*normal,2);
tangent = c.optics.nWater*(incident-incidence.*normal);
q = vecnorm(tangent,2,2);
tir = q>=c.optics.nAir;
T = rp_fresnel(c.optics.nWater,c.optics.nAir,[q,zeros(size(q))]);
T(incidence<=0) = 0;
kz = sqrt(max(0,c.optics.nAir^2-q.^2));
direction = (tangent+kz.*normal)/c.optics.nAir;
end

function d = scatterHG(d,g)
% Same sampler as W_RIS/wris_mc (yimiao Hen_Green for the default g = 0.924).
n = size(d,1);
if n==0, return; end
if g==0.924
    mu = real(cos(Hen_Green(n,2).'));
else
    u = rand(n,1);
    if abs(g)<1e-12
        mu = 2*u-1;
    else
        mu = (1+g^2-((1-g^2)./(1-g+2*g*u)).^2)/(2*g);
    end
end
mu = max(-1,min(1,mu));
azimuth = 2*pi*rand(n,1);
helper = repmat([0 0 1],n,1);
polar = abs(d(:,3))>0.9;
helper(polar,:) = repmat([1 0 0],nnz(polar),1);
e1 = cross(helper,d,2);
e1 = e1./vecnorm(e1,2,2);
e2 = cross(d,e1,2);
d = mu.*d+sqrt(max(0,1-mu.^2)).*(cos(azimuth).*e1+sin(azimuth).*e2);
d = d./vecnorm(d,2,2);
end
