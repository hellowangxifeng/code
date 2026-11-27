function result=rp_fixed_refine_unmodulated(c,d,previous,samples)
%RP_FIXED_REFINE_UNMODULATED Smooth empirical incoherent scatter consistently.
% Scott bandwidth is set from packet effective count, never from BER ranking.
weights=samples.powerW; xy=samples.xy; total=sum(weights); normalized=weights/total;
center=sum(xy.*normalized,1); delta=xy-center;
covariance=(delta.*normalized).'*delta;
effectiveCount=total^2/sum(weights.^2);
bandwidth=sqrt(trace(covariance)/2)*effectiveCount^(-1/6);
if isfield(samples,'pilotBandwidthM'),bandwidth=samples.pilotBandwidthM;end
target=previous.receiverCenterXY; kinds=cell(1,3);
for k=1:3
    use=true(size(weights)); factor=1;
    if k==2, use=samples.blockId<=samples.blocks/2; factor=2; end
    if k==3, use=samples.blockId>samples.blocks/2; factor=2; end
    wave=scatterWave(xy(use,:),factor*weights(use),target,bandwidth);
    channel=rp_aperture_channel(c,d.geometry,wave,target,rp_link_samples(c));
    kinds{k}=struct('wave',wave,'channel',channel);
end
powerW=10.^((0:.5:100)/10)/1000;
first=rp_fixed_receiver_metrics(c,kinds{2}.channel,powerW);
second=rp_fixed_receiver_metrics(c,kinds{3}.channel,powerW);
errors=[max(abs(first.aber-second.aber)),max(abs(first.op-second.op))];
fprintf('Independent MC halves: ABER %.8g, OP %.8g; Scott bandwidth %.6g m\n',errors(1),errors(2),bandwidth);
assert(errors(1)<.015 && errors(2)<.03,'rp:FixedScatterConvergence','Independent scatter MC halves disagree.');
wave=kinds{1}.wave; result=previous;
result.channelLower=kinds{1}.channel;
wave.irradianceWm2=wave.irradianceWm2+previous.tailIrradianceEnvelopeWm2;
result.channelUpper=rp_aperture_channel(c,d.geometry,wave,target,rp_link_samples(c));
result.wave=wave; result.scatterPhotons=samples.photons; result.scatterEffectivePackets=effectiveCount;
result.scatterPacketsInWindow=numel(weights); result.scatterPowerInWindowW=total;
result.scatterBandwidthM=bandwidth; result.scatterIndependentABERDifference=errors(1);
result.scatterIndependentOPDifference=errors(2);
result.scatterEstimator=sprintf('Weighted pilot-Scott Gaussian kernel; %d independent photon batches; two independent half-runs',samples.blocks);
% A bandwidth sensitivity check quantifies smoothing, rather than displaying
% empirical empty histogram pixels as a physical zero-light floor.
sensitivity=zeros(3,3);
for k=1:3
    scale=[.75 1 1.25]; w=scatterWave(xy,weights,target,scale(k)*bandwidth);
    h=rp_aperture_channel(c,d.geometry,w,target,rp_link_samples(c));
    m=rp_fixed_receiver_metrics(c,h,powerW); reference=rp_fixed_receiver_metrics(c,kinds{1}.channel,powerW);
    sensitivity(k,:)=[scale(k),max(abs(m.aber-reference.aber)),max(abs(m.op-reference.op))];
end
result.scatterBandwidthSensitivity=sensitivity;
assert(max(sensitivity(:,2))<.015 && max(sensitivity(:,3))<.03, ...
    'rp:FixedScatterBandwidth','Off-axis scatter smoothing is not converged.');
save(fullfile(c.output,'fixed_receiver','unmodulated_refined.mat'),'result','-v7.3');
end

function w=scatterWave(xy,power,target,bandwidth)
half=.65; n=192;
xe=linspace(target(1)-half,target(1)+half,n+1); ye=linspace(target(2)-half,target(2)+half,n+1);
dx=xe(2)-xe(1); dy=ye(2)-ye(1);
ix=discretize(xy(:,1),xe); iy=discretize(xy(:,2),ye); use=isfinite(ix) & isfinite(iy);
raw=accumarray([iy(use),ix(use)],power(use),[n n],@sum,0)/(dx*dy);
t=(-ceil(4*bandwidth/dx):ceil(4*bandwidth/dx))*dx;
kernel=exp(-t.^2/(2*bandwidth^2)); kernel=kernel/sum(kernel);
I=conv2(kernel,kernel,raw,'same');
w=struct('x',(xe(1:end-1)+xe(2:end))/2,'y',(ye(1:end-1)+ye(2:end))/2, ...
    'irradianceWm2',I,'cellAreaM2',dx*dy);
end
