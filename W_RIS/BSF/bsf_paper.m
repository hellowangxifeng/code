function result = bsf_paper(c, g, x, y, options)
%BSF_PAPER Equation (3-18), separated analytically into Ib and Is.
% k1=1.05,k2=-0.07973 are the UNFITTED constants in fig3-10.
% The finite waist is matched to W_RIS; no air-water leg is introduced.
arguments
    c (1,1) struct
    g (1,1) struct
    x {mustBeNumeric, mustBeReal, mustBeFinite}
    y {mustBeNumeric, mustBeReal, mustBeFinite}
    options.FrequencyNodes (1,1) {mustBeInteger, mustBePositive, mustBeReal} = 2049
    options.K1 (1,1) {mustBePositive, mustBeReal} = 1.05
    options.K2 (1,1) {mustBeNegative, mustBeReal} = -0.07973
end
assert(isequal(size(x),size(y)), 'bsf:Shape', 'x and y must have matching sizes.');
assert(mod(options.FrequencyNodes,2)==1 && options.FrequencyNodes>=3, ...
    'bsf:FrequencyNodes', 'Composite Simpson integration needs an odd number >=3.');
q = [x(:),y(:),zeros(numel(x),1)]*g.RRis.' + g.pRis.';
r = (q-g.pTx.')*g.RTx;
z = r(:,3);
rho = vecnorm(r(:,1:2),2,2);
mu = max((q-g.pTx.')*g.normal,0)./vecnorm(r,2,2);
valid = z>0;
width = sqrt(c.source.waist^2+(c.source.divergence*z).^2);
extinction = c.water.absorption+c.water.scattering;
Ib = zeros(size(z));
Ib(valid) = c.source.power*2./(pi*width(valid).^2).* ...
    exp(-2*rho(valid).^2./width(valid).^2-extinction*z(valid)).*mu(valid);
Is = zeros(size(z));
kmax = sqrt(8*log(1e14))/min(width(valid));
k = linspace(0,kmax,options.FrequencyNodes);
dk = k(2)-k(1);
simpson = ones(size(k));
simpson(2:2:end-1) = 4;
simpson(3:2:end-1) = 2;
simpson = simpson*dk/3;
indices = find(valid);
for start = 1:256:numel(indices)
    ids = indices(start:min(start+255,numel(indices)));
    argument = options.K2*z(ids).*k;
    % expm1(t)/t -> 1 exactly at t=0, avoiding the original 0/0.
    ratio = ones(size(argument));
    ratio(:,2:end) = expm1(argument(:,2:end))./argument(:,2:end);
    H = c.water.scattering*options.K1*z(ids).*ratio;
    spectralScatter = exp(-width(ids).^2.*k.^2/8).*expm1(H);
    Is(ids) = c.source.power/(2*pi)*exp(-extinction*z(ids)).*mu(ids).* ...
        sum(besselj(0,rho(ids).*k).*spectralScatter.*(k.*simpson),2);
end
result = struct('Ib',reshape(Ib,size(x)),'Is',reshape(Is,size(x)), ...
    'IBSF',reshape(Ib+Is,size(x)), 'kmax',kmax, ...
    'frequencyNodes',options.FrequencyNodes,'k1',options.K1,'k2',options.K2, ...
    'minimumScatter',min(Is));
end
