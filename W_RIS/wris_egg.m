function [h,info] = wris_egg(kind,n,seed)
%WRIS_EGG Mean-one EGG irradiance gains from supplied image Table 2-1.
% [h,info]=wris_egg('weak'|'strong'|'none',n,seed)
% n denotes channel realizations/time snapshots, NOT individual photons.
% Uses randg (Statistics and Machine Learning Toolbox in MATLAB R2024b).
% X ~ w*Exp(scale=lambda)+(1-w)*GG(a,b,c), H=X/E[X].
% A snapshot shares one h over its entire photon cloud / received map.
% The supplied table specifies a marginal law, not a spatial covariance,
% temporal correlation, path-length scaling, or stochastic ray deflections.

arguments
    kind = 'weak'
    n (1,1) {mustBeNumeric,mustBeReal,mustBeFinite,mustBeNonnegative,mustBeInteger} = 1
    seed {mustBeSeed} = []
end
kind = validatestring(kind,{'weak','strong','none'},mfilename,'kind');
if ~isempty(seed)
    oldRng = rng;
    cleanup = onCleanup(@() rng(oldRng));
    rng(seed,'twister');
end
info.kind = kind;
info.sampleCount = n;
info.seed = seed;
info.application = 'One gain per entire channel snapshot, shared by all photons.';
info.momentOrders = 0:4;
if strcmp(kind,'none')
    h = ones(n,1);
    info.parameters = struct('omega',0,'lambda_e',NaN,'a_g',NaN,'b_g',NaN,'c_g',NaN);
    info.rawMean = 1;
    info.rawMoments = ones(1,5);
    info.normalizedMoments = ones(1,5);
    info.theoreticalScintillationIndex = 0;
    info.tabulatedScintillationIndex = 0;
    info.bubbleFlow_Lmin = 0;
    info.temperatureGradient_CperCm = 0;
    return;
end
if strcmp(kind,'weak')
    p = [0.2130 0.3291 1.4299 1.1817 17.1984];
    tabSI = 0.1484; flow = 2.4; gradient = 0.05;
else
    p = [0.4539 0.2744 0.3008 1.7053 54.1422];
    tabSI = 0.4769; flow = 4.7; gradient = 0.1;
end
w=p(1); lambda=p(2); a=p(3); b=p(4); c=p(5);
k = info.momentOrders;
rawMoments = w*exp(k*log(lambda)+gammaln(k+1)) + ...
    (1-w)*exp(k*log(b)+gammaln(a+k/c)-gammaln(a));
mu = rawMoments(2);
isExp = rand(n,1)<w;
x = zeros(n,1);
ne = nnz(isExp);
x(isExp) = -lambda*log(max(rand(ne,1),realmin));
ng = n-ne;
if ng > 0
    x(~isExp) = b*randg(a,ng,1).^(1/c);
end
h = x/mu;
info.parameters = struct('omega',w,'lambda_e',lambda,'a_g',a,'b_g',b,'c_g',c);
info.rawMean = mu;
info.rawMoments = rawMoments;
info.normalizedMoments = rawMoments./mu.^k;
info.theoreticalScintillationIndex = rawMoments(3)/mu^2-1;
info.tabulatedScintillationIndex = tabSI;
info.bubbleFlow_Lmin = flow;
info.temperatureGradient_CperCm = gradient;
info.source = 'Supplied image-1 Eqs. (2-10),(2-11); image-2 Eqs. (2-12),(2-13), Table 2-1.';
end

function mustBeSeed(seed)
if isempty(seed)
    mustBeNumeric(seed);
else
    validateattributes(seed,{'numeric'}, ...
        {'scalar','real','finite','nonnegative','integer','<=',2^32-1});
end
end
