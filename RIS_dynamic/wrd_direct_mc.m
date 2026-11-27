function out = wrd_direct_mc(c,g,sea,seed)
%WRD_DIRECT_MC Implicit-capture HG transport with a rough water-air boundary.
old = rng;
restore = onCleanup(@() rng(old));
rng(seed,'twister');
photons = c.mc.photons;
absorption = c.water.absorption;
scattering = c.water.scattering;
extinction = absorption+scattering;
ledger = struct('incident',0,'fresnelReflected',0,'tir',0,'transmitted',0, ...
    'absorbed',0,'backEscape',0,'truncated',0,'rouletteDelta',0,'waveBlocked',0);
chunks = cell(ceil(photons/c.mc.chunkSize),1);
for first = 1:c.mc.chunkSize:photons
    n = min(c.mc.chunkSize,photons-first+1);
    uv = c.source.waist/2*randn(n,2);
    slope = c.source.divergence/2*randn(n,2);
    launch = [slope,ones(n,1)];
    direction = (launch./vecnorm(launch,2,2))*g.RTx.';
    position = [uv,zeros(n,1)]*g.RTx.'+g.pTx.';
    weights = ones(n,1);
    path = zeros(n,1);
    collisions = zeros(n,1);
    active = true(n,1);
    airPosition = zeros(n,3);
    airDirection = zeros(n,3);
    airWeight = zeros(n,1);
    % Test launch medium only when the height bound cannot certify it.
    nearSurface = position(:,3)>=-sea.amplitudeBound;
    if any(nearSurface)
        assert(all(position(nearSurface,3)<sea.evaluate(position(nearSurface,1:2))), ...
            'wrd:SourceMedium','The TX launch plane must start under water.');
    end
    while any(active)
        ids = find(active);
        p = position(ids,:);
        d = direction(ids,:);
        free = inf(size(ids));
        if extinction>0
            free = -log(max(rand(size(ids)),realmin))/extinction;
        end
        remaining = c.mc.maxPath-path(ids);
        % Use the same source-plane back-escape boundary as W_RIS.
        axial = d*g.dTx;
        back = ((g.pTx.'-p)*g.dTx)./axial;
        back(axial>=-1e-12 | back<1e-10) = Inf;
        limit = min([free,remaining,back],[],2);
        surface = wrd_crossing(p,d,limit,sea,c.surface,1);
        hit = isfinite(surface);
        reached = ids(hit);
        if any(hit)
            exit = p(hit,:)+d(hit,:).*surface(hit);
            [~,gradient] = sea.evaluate(exit(:,1:2));
            normal = [-gradient,ones(nnz(hit),1)];
            normal = normal./vecnorm(normal,2,2);
            [airDirection(reached,:),transmission,tir] = refract(c,d(hit,:),normal);
            airPosition(reached,:) = exit;
            airWeight(reached) = weights(reached).*transmission;
            ledger.incident = ledger.incident+sum(weights(reached));
            ledger.tir = ledger.tir+sum(weights(reached(tir)));
            ledger.transmitted = ledger.transmitted+sum(airWeight(reached));
            ledger.fresnelReflected = ledger.fresnelReflected+ ...
                sum(weights(reached(~tir)).*(1-transmission(~tir)));
        end
        escaped = ~hit & back<=free & back<=remaining;
        capped = ~hit & ~escaped & remaining<=free;
        ledger.backEscape = ledger.backEscape+sum(weights(ids(escaped)));
        ledger.truncated = ledger.truncated+sum(weights(ids(capped)));
        active(ids(hit | escaped | capped)) = false;
        move = ~hit & ~escaped & ~capped;
        collide = ids(move);
        if isempty(collide), continue; end
        position(collide,:) = p(move,:)+d(move,:).*free(move);
        path(collide) = path(collide)+free(move);
        collisions(collide) = collisions(collide)+1;
        ledger.absorbed = ledger.absorbed+sum(weights(collide))*absorption/extinction;
        weights(collide) = weights(collide)*scattering/extinction;
        low = collide(weights(collide)>0 & weights(collide)<c.mc.rouletteThreshold);
        before = weights(low);
        survives = rand(size(low))<c.mc.rouletteSurvival;
        weights(low) = before.*survives/c.mc.rouletteSurvival;
        ledger.rouletteDelta = ledger.rouletteDelta+sum(weights(low)-before);
        active(collide(weights(collide)==0)) = false;
        killed = collide(active(collide) & collisions(collide)>=c.mc.maxScatter);
        ledger.truncated = ledger.truncated+sum(weights(killed));
        active(killed) = false;
        scatter = collide(active(collide));
        direction(scatter,:) = scatterHG(direction(scatter,:),c.water.anisotropy);
    end
    transmitted = airWeight>0;
    [arrivals,blocked] = wrd_receive(c,airPosition(transmitted,:), ...
        airDirection(transmitted,:),airWeight(transmitted),collisions(transmitted)==0,sea);
    ledger.waveBlocked = ledger.waveBlocked+blocked;
    chunks{ceil(first/c.mc.chunkSize)} = arrivals;
end
out = wrd_score(c,vertcat(chunks{:}));
scale = c.source.power/photons;
out.transmittedPowerW = ledger.transmitted*scale;
out.tirPowerW = ledger.tir*scale;
out.incidentPowerW = ledger.incident*scale;
out.waveBlockedPowerW = ledger.waveBlocked*scale;
out.energyResidual = (ledger.transmitted+ledger.fresnelReflected+ledger.tir+ ...
    ledger.absorbed+ledger.backEscape+ledger.truncated-photons-ledger.rouletteDelta)/photons;
out.ledger = ledger;
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
n = size(d,1);
if n==0, return; end
if g==0.924
    mu = cos(Hen_Green(n,2).');
else
    u = rand(n,1);
    if abs(g)<1e-12
        mu = 2*u-1;
    else
        mu = (1+g^2-((1-g^2)./(1-g+2*g*u)).^2)/(2*g);
    end
end
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
