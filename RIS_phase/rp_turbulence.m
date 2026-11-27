function t=rp_turbulence(c,order)
%RP_TURBULENCE Mean-one ocean EGG times independent atmospheric Gamma-Gamma.
% GG = Gamma(alpha,1/alpha)*Gamma(beta,1/beta); shapes are explicit scenarios.
if nargin<2, order=32; end
p=parameters(c.study.turbulence);
[x,we]=laguerre(order,1); [y,wg]=laguerre(order,p(3));
if strcmp(c.study.turbulence,'none')
    ocean=1; wo=1;
else
    ocean=[p(2)*x;p(4)*y.^(1/p(5))]/p(6);
    wo=[p(1)*we;(1-p(1))*wg];
end
kind='none'; if isfield(c.study,'atmosphere'), kind=c.study.atmosphere; end
switch kind
    case 'none', atmosphere=1; wa=1; shapes=[Inf Inf];
    case 'weak', shapes=c.study.ggWeak;
    case 'strong', shapes=c.study.ggStrong;
    otherwise, error('rp:Atmosphere','Unknown atmospheric scenario.');
end
if ~strcmp(kind,'none')
    [a,pa]=laguerre(order,shapes(1)); [b,pb]=laguerre(order,shapes(2));
    atmosphere=reshape((a/shapes(1))*(b/shapes(2)).',[],1);
    wa=reshape(pa*pb.',[],1);
end
values=ocean*atmosphere.'; weights=wo*wa.';
t=struct('values',values(:),'weights',weights(:),'oceanParameters',p, ...
    'atmosphereValues',atmosphere,'atmosphereWeights',wa,'ggShapes',shapes, ...
    'ocean',c.study.turbulence,'atmosphere',kind,'order',order);
end

function p=parameters(kind)
switch kind
    case 'none', p=[0 1 1 1 1 1];
    case 'weak', p=[.2130 .3291 1.4299 1.1817 17.1984 0];
    case 'strong', p=[.4539 .2744 .3008 1.7053 54.1422 0];
    otherwise, error('rp:Turbulence','Unknown ocean EGG scenario.');
end
if p(6)==0
    p(6)=p(1)*p(2)+(1-p(1))*p(4)*exp(gammaln(p(3)+1/p(5))-gammaln(p(3)));
end
end

function [x,w]=laguerre(n,a)
i=(1:n).'; off=sqrt((1:n-1).'.*((1:n-1).'+a-1));
[V,D]=eig(diag(2*i+a-2)+diag(off,1)+diag(off,-1));
[x,idx]=sort(diag(D)); w=V(1,idx).'.^2; w=w/sum(w);
end
