function s=rp_link_samples(c,n,seed)
%RP_LINK_SAMPLES Independent randomized Sobol samples; untruncated Gaussian tails.
if nargin<2, n=c.study.samples; end
if nargin<3, seed=c.study.seed; end
old=rng; cleanup=onCleanup(@() rng(old)); rng(seed,'twister');
q=scramble(sobolset(4,'Skip',1024),'MatousekAffineOwen');
u=net(q,n); z=sqrt(2)*erfinv(2*max(eps,min(1-eps,u))-1);
s=struct('positionNormal',z(:,1:2),'angleNormal',z(:,3:4),'n',n,'seed',seed);
end
