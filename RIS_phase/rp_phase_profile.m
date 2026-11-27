function f = rp_phase_profile(c,p,m,delta)
%RP_PHASE_PROFILE One scalar phase and its analytic, curl-free gradient.
% m may be one unit or one unit index for every coordinate row.
arguments
    c (1,1) struct
    p (1,1) struct
    m (:,1) double {mustBeInteger,mustBePositive}
    delta (:,2) double {mustBeReal,mustBeFinite}
end
if isscalar(m), m=repmat(m,size(delta,1),1); end
assert(numel(m)==size(delta,1),'rp:ProfileSize','Provide one unit index per point.');
b=p.gradientOverK0(m,:);
h11=reshape(p.phaseCurvature(1,1,:),[],1);
h12=reshape(p.phaseCurvature(1,2,:),[],1);
h22=reshape(p.phaseCurvature(2,2,:),[],1);
model='quadratic';
if isfield(p,'phaseModel'), model=p.phaseModel; end
if strcmp(model,'tanh')
    assert(all(abs(h12)<1e-14),'rp:Integrability','The separable tanh phase requires Hxy=0.');
    B=c.optics.nAir*sind(c.device.maxAngleDeg);
    H=[h11(m),h22(m)];
    assert(all(abs(b)<B,'all'),'rp:TanhDomain','The tanh reference gradient must satisfy abs(b)<B.');
    a=H./(B*(1-(b/B).^2));
    t=a.*delta;
    z=atanh(b/B)+t;
    gradient=B*tanh(z);
    % Stable log(cosh(c+t)/cosh(c)); avoid overflow and small-H cancellation.
    difference=logCosh(z)-logCosh(atanh(b/B));
    small=abs(t)<0.5;
    difference(small)=log1p(2*sinh(t(small)/2).^2+(b(small)/B).*sinh(t(small)));
    value=B*difference./a;
    tiny=abs(t)<1e-6;
    series=b.*delta+H.*delta.^2/2-a.*H.*(b/B).*delta.^3/3;
    value(tiny)=series(tiny);
    localH=B*a.*(1-tanh(z).^2);
    phase=sum(value,2);
else
    gradient=b+[h11(m).*delta(:,1)+h12(m).*delta(:,2), ...
        h12(m).*delta(:,1)+h22(m).*delta(:,2)];
    phase=sum(b.*delta,2)+(h11(m).*delta(:,1).^2+ ...
        2*h12(m).*prod(delta,2)+h22(m).*delta(:,2).^2)/2;
    localH=[h11(m),h22(m)];
end
f=struct('phaseOverK0',phase,'gradientOverK0',gradient,'diagonalH',localH);
end

function y=logCosh(x)
y=abs(x)+log1p(exp(-2*abs(x)))-log(2);
end
