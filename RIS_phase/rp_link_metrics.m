function m=rp_link_metrics(c,h,powerW,order)
%RP_LINK_METRICS Joint EGG/GG fading, thermal OOK and finite-aperture jitter.
if nargin<4, order=32; end
response=rp_fading_response(c,order); gain=h.gain(:); powerW=powerW(:).';
P=gain*powerW; logP=log10(max(P,realmin));
R=c.study.responsivityAW; noise=c.study.noiseStdA;
if response.deterministic
    ber=.5*erfc(R*P/(2*sqrt(2)*noise)); outage=double(P<response.threshold);
else
    logP=max(response.logPower(1),min(response.logPower(end),logP));
    ber=interp1(response.logPower,response.ber,logP,'pchip');
    outage=interp1(response.logPower,response.cdf,logP,'pchip');
    ber(P==0)=.5; outage(P==0)=1;
end
m=struct('powerW',powerW,'meanReceivedPowerW',mean(gain)*powerW, ...
    'op',mean(outage,1),'aber',mean(ber,1),'fovOutageExact',h.fovOutageExact, ...
    'fovOutageSample',h.fovOutageSample,'order',order, ...
    'thresholdReceivedPowerW',response.threshold,'fading','independent mean-one EGG x Gamma-Gamma');
end
