function tx = wris_tx_motion(c,t)
%WRIS_TX_MOTION Weng et al. (2022), Eqs (5)-(7): planar dead reckoning.
% The paper supplies measured u,v,r, not a stochastic law. OU inputs here are
% explicitly uncalibrated scenario assumptions. No PF/PAT controller assumed.
ref=wris_reference(c); n=numel(t);
stream=RandStream('mt19937ar','Seed',c.tx.seed);
uvr=zeros(3,n); uvr(:,1)=c.tx.meanUVR+c.tx.stdUVR.*randn(stream,3,1);
p=repmat(ref.pTx,1,n); yaw=zeros(1,n); yaw(1)=c.tx.initialYaw;
for k=2:n
    dt=t(k)-t(k-1);
    q=exp(-dt./c.tx.correlationTime);
    uvr(:,k)=c.tx.meanUVR+q.*(uvr(:,k-1)-c.tx.meanUVR)+ ...
        c.tx.stdUVR.*sqrt(1-q.^2).*randn(stream,3,1);
    psi=yaw(k-1); u=uvr(1,k-1); v=uvr(2,k-1);
    p(1:2,k)=p(1:2,k-1)+dt*[u*cos(psi)-v*sin(psi);u*sin(psi)+v*cos(psi)];
    yaw(k)=psi+dt*uvr(3,k-1);
end
rotation=zeros(3,3,n); d=zeros(3,n);
for k=1:n
    if ~strcmp(c.tx.pointing,'body-fixed'), error('wris:Pointing','Only body-fixed pointing is supported.'); end
    R=RIS_rotation([0;0;yaw(k)-yaw(1)]);
    rotation(:,:,k)=R*ref.RTx; d(:,k)=rotation(:,3,k);
end
tx=struct('t',t,'position',p,'yaw',yaw,'uvr',uvr,'rotation',rotation,'direction',d, ...
    'model','111.pdf Eqs (5)-(7), constant depth, assumed OU velocity inputs');
end
