function [out, packets] = bsf_scatter_transport(c, g, options)
%BSF_SCATTER_TRANSPORT Unbiased evaluation of the >=1-collision HG path integral.
% Force the first scattering before the first scoring boundary; subsequently
% use scattering-only free paths and continuous absorption. No MC-map fitting.
% This computes Is; it is NOT an independent closed-form approximation to Is.
arguments
    c (1,1) struct
    g (1,1) struct
    options.Photons (1,1) {mustBeNumeric, mustBeInteger, mustBeReal, mustBePositive} = 2e6
    options.Seed (1,1) {mustBeNumeric, mustBeInteger, mustBeReal, mustBeNonnegative} = 73101
    options.SavePackets (1,1) logical = false
    options.Sampling (1,1) string {mustBeMember(options.Sampling,["pseudo","sobol"])} = "pseudo"
end
old = rng;
guard = onCleanup(@() rng(old));
rng(options.Seed,'twister');
N = options.Photons;
assert(N>1, 'bsf:Photons', 'At least two samples are needed for a standard error.');
assert(c.water.absorption>=0 && c.water.scattering>=0, 'bsf:Water', 'a and b must be nonnegative.');
nx = c.grid.bins(1);
ny = c.grid.bins(2);
xe = linspace(-c.ris.Lx/2,c.ris.Lx/2,nx+1);
ye = linspace(-c.ris.Ly/2,c.ris.Ly/2,ny+1);
area = (xe(2)-xe(1))*(ye(2)-ye(1));
h = zeros(nx,ny);
h2 = h;
b = c.water.scattering;
a = c.water.absorption;
numChunks = ceil(N/c.mc.chunkSize);
packetCells = cell(numChunks,1);
powerSum = 0;
powerSum2 = 0;
truncatedSum = 0;
acceptedCount = 0;
useSobol = options.Sampling=="sobol";
if useSobol
    assert(exist('sobolset','file')==2, 'bsf:Sobol', 'Sobol sampling requires Statistics and Machine Learning Toolbox.');
    sobol = scramble(sobolset(32),'MatousekAffineOwen');
end
for first = 1:c.mc.chunkSize:N
    if b==0
        break
    end
    chunk = ceil(first/c.mc.chunkSize);
    n = min(c.mc.chunkSize,N-first+1);
    uniforms = [];
    if useSobol
        uniforms = sobol(first:first+n-1,:);
        gaussian = sqrt(2)*erfinv(2*uniforms(:,1:4)-1);
        uv = c.source.waist/2*gaussian(:,1:2);
        slope = c.source.divergence/2*gaussian(:,3:4);
    else
        uv = c.source.waist/2*randn(n,2);
        slope = c.source.divergence/2*randn(n,2);
    end
    d = [slope,ones(n,1)];
    d = (d./vecnorm(d,2,2))*g.RTx.';
    p = [uv,zeros(n,1)]*g.RTx.' + g.pTx.';
    [~,planeDistance] = wris_intersect(p,d,g.pRis,g.RRis,[c.ris.Lx,c.ris.Ly]);
    planeDistance(d*g.normal<=1e-12) = inf;
    limit = min(planeDistance,c.mc.maxPath);
    forceProbability = -expm1(-b*limit);
    pathLength = -log1p(-getUniforms(uniforms,(1:n).',5,1).*forceProbability)/b;
    p = p + d.*pathLength;
    w = forceProbability.*exp(-a*pathLength);
    d = scatterDirection(d,c.water.anisotropy,getUniforms(uniforms,(1:n).',6,2));
    active = true(n,1);
    order = ones(n,1);
    chunkPackets = cell(0,1);
    packetIndex = 0;
    iteration = 0;
    while any(active)
        iteration = iteration+1;
        ids = find(active);
        pp = p(ids,:);
        dd = d(ids,:);
        dimension = 8+3*(iteration-1);
        free = -log(max(getUniforms(uniforms,ids,dimension,1),realmin))/b;
        remaining = c.mc.maxPath-pathLength(ids);
        travel = min(free,remaining);
        [local,dist] = wris_intersect(pp,dd,g.pRis,g.RRis,[c.ris.Lx,c.ris.Ly]);
        dist(dd*g.normal<=1e-12) = inf;
        denom = dd*g.dTx;
        back = ((g.pTx.'-pp)*g.dTx)./denom;
        back(denom>=-1e-12 | back<1e-10) = inf;
        reach = dist<=travel & dist<=back & isfinite(dist);
        escape = back<=travel & back<dist & isfinite(back);
        hits = local(reach,1:2);
        inside = abs(hits(:,1))<=c.ris.Lx/2 & abs(hits(:,2))<=c.ris.Ly/2;
        hitWeights = w(ids(reach)).*exp(-a*dist(reach));
        hitWeights = hitWeights(inside);
        hitDistance = dist(reach);
        hitDistance = hitDistance(inside);
        hits = hits(inside,:);
        if ~isempty(hitWeights)
            ix = discretize(hits(:,1),xe);
            iy = discretize(hits(:,2),ye);
            h = h + accumarray([ix,iy],hitWeights,[nx,ny]);
            h2 = h2 + accumarray([ix,iy],hitWeights.^2,[nx,ny]);
        end
        powerSum = powerSum + sum(hitWeights);
        powerSum2 = powerSum2 + sum(hitWeights.^2);
        acceptedCount = acceptedCount + numel(hitWeights);
        if options.SavePackets && ~isempty(hitWeights)
            hitIds = ids(reach);
            hitIds = hitIds(inside);
            packetIndex = packetIndex+1;
            chunkPackets{packetIndex,1} = [hits,zeros(size(hits,1),1),d(hitIds,:), ...
                hitWeights*c.source.power/N,pathLength(hitIds)+hitDistance];
        end
        active(ids(reach|escape)) = false;
        moveIds = ids(~reach & ~escape);
        move = ~reach & ~escape;
        if isempty(moveIds)
            continue
        end
        p(moveIds,:) = pp(move,:) + dd(move,:).*travel(move);
        pathLength(moveIds) = pathLength(moveIds) + travel(move);
        w(moveIds) = w(moveIds).*exp(-a*travel(move));
        cutoff = free(move)>=remaining(move);
        cutIds = moveIds(cutoff);
        truncatedSum = truncatedSum + sum(w(cutIds));
        active(cutIds) = false;
        collisionIds = moveIds(~cutoff);
        order(collisionIds) = order(collisionIds)+1;
        cutIds = collisionIds(order(collisionIds)>=c.mc.maxScatter);
        truncatedSum = truncatedSum + sum(w(cutIds));
        active(cutIds) = false;
        collisionIds = collisionIds(active(collisionIds));
        d(collisionIds,:) = scatterDirection(d(collisionIds,:),c.water.anisotropy, ...
            getUniforms(uniforms,collisionIds,dimension+1,2));
    end
    if options.SavePackets && ~isempty(chunkPackets)
        packetCells{chunk} = vertcat(chunkPackets{:});
    end
end
scale = c.source.power/N;
variance = max(0,h2-h.^2/N)/(N*(N-1));
out = struct('Is',h*scale/area,'standardError',sqrt(variance)*c.source.power/area, ...
    'x',(xe(1:end-1)+xe(2:end))/2,'y',(ye(1:end-1)+ye(2:end))/2, ...
    'power',powerSum*scale,'powerSE',c.source.power*sqrt(max(0,powerSum2-powerSum^2/N)/(N*(N-1))), ...
    'truncatedPower',truncatedSum*scale,'acceptedCount',acceptedCount, ...
    'photons',N,'seed',options.Seed,'cellArea',area, ...
    'sampling',options.Sampling, ...
    'method','forced-first-scatter, continuous absorption, exact HG');
if useSobol
    out.iidSEProxy = out.standardError;
    out.standardError = []; % Must use BETWEEN independent scrambles, never iid variance/N.
    out.powerSE = NaN;
end
packets = struct('localXY',[],'pos',[],'dir',[],'powerW',[],'pathLengthM',[]);
if options.SavePackets
    rows = vertcat(packetCells{:});
    if ~isempty(rows)
        packets.localXY = rows(:,1:2);
        packets.pos = rows(:,1:3)*g.RRis.' + g.pRis.';
        packets.dir = rows(:,4:6);
        packets.powerW = rows(:,7);
        packets.pathLengthM = rows(:,8);
    end
end
end

function uniform = getUniforms(values,ids,dimension,count)
if ~isempty(values) && dimension+count-1<=size(values,2)
    uniform = values(ids,dimension:dimension+count-1);
else
    uniform = rand(numel(ids),count);
end
end

function d = scatterDirection(d,g,uniform)
n = size(d,1);
u = uniform(:,1);
if abs(g)<1e-12
    mu = 2*u-1;
else
    mu = (1+g^2-((1-g^2)./(1-g+2*g*u)).^2)/(2*g);
end
mu = max(-1,min(1,mu));
phi = 2*pi*uniform(:,2);
helper = repmat([0,0,1],n,1);
polar = abs(d(:,3))>0.9;
helper(polar,:) = repmat([1,0,0],nnz(polar),1);
e1 = cross(helper,d,2);
e1 = e1./vecnorm(e1,2,2);
e2 = cross(d,e1,2);
d = mu.*d + sqrt(1-mu.^2).*(cos(phi).*e1+sin(phi).*e2);
d = d./vecnorm(d,2,2);
end
