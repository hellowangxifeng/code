function response=rp_fading_response(c,order)
%RP_FADING_RESPONSE Tail-resolving log-probability convolution of ocean EGG*GG.
% Exact CDF bin masses avoid the deep-fade bias of low-order Laguerre BER
% integration. order controls log-grid refinement (32:0.0025, 64:0.00125).
if nargin<2, order=32; end
persistent keys tables
if isempty(keys), keys={}; tables={}; end
kind='none'; if isfield(c.study,'atmosphere'), kind=c.study.atmosphere; end
key=jsonencode({c.study.turbulence,kind,order,c.study.responsivityAW, ...
    c.study.noiseStdA,c.study.snrThreshold});
if ~strcmp(kind,'none'), key=[key,jsonencode({c.study.ggWeak,c.study.ggStrong})]; end
idx=find(strcmp(keys,key),1);
if ~isempty(idx), response=tables{idx}; return; end
parameters=rp_turbulence(c,2); p=parameters.oceanParameters;
logPower=linspace(-13,3,4097).'; power=10.^logPower;
R=c.study.responsivityAW; noise=c.study.noiseStdA;
threshold=sqrt(c.study.snrThreshold)*noise/R;
if strcmp(c.study.turbulence,'none') && strcmp(kind,'none')
    response=struct('logPower',logPower,'ber',.5*erfc(R*power/(2*sqrt(2)*noise)), ...
        'cdf',double(power<threshold),'deterministic',true,'threshold',threshold, ...
        'logStep',0,'distributionMass',1,'distributionMean',1);
else
    step=.08/order; edges=-45:step:12; centers=(edges(1:end-1)+edges(2:end))/2;
    mass=1; origin=0;
    if ~strcmp(c.study.turbulence,'none')
        h=exp(edges);
        F=p(1)*(-expm1(-p(6)*h/p(2)))+ ...
            (1-p(1))*gammainc((p(6)*h/p(4)).^p(5),p(3),'lower');
        [part,start]=trimMass(max(0,diff(F)),centers);
        mass=part; origin=start;
    end
    if ~strcmp(kind,'none')
        for shape=parameters.ggShapes
            F=gammainc(shape*exp(edges),shape,'lower');
            [part,start]=trimMass(max(0,diff(F)),centers);
            [mass,origin]=convolveMass(mass,origin,part,start);
        end
    end
    logH=origin+(0:numel(mass)-1)*step;
    mass=max(0,mass); mass(mass<1e-18)=0; probability=sum(mass); mass=mass/probability;
    meanH=sum(exp(logH).*mass);
    distributionCDF=cumsum(mass)-mass/2;
    cdf=interp1(logH,distributionCDF,log(threshold./power),'pchip','extrap');
    cdf=max(0,min(1,cdf));
    keep=mass>0; values=exp(logH(keep)); weights=mass(keep).'; ber=zeros(size(power));
    for first=1:256:numel(values)
        ids=first:min(first+255,numel(values));
        ber=ber+(.5*erfc((R/(2*sqrt(2)*noise))*power*values(ids)))*weights(ids);
    end
    response=struct('logPower',logPower,'ber',ber,'cdf',cdf,'deterministic',false, ...
        'threshold',threshold,'logStep',step,'distributionMass',probability, ...
        'distributionMean',meanH);
end
keys{end+1}=key; tables{end+1}=response;
end

function [mass,origin]=trimMass(mass,nodes)
keep=find(mass>1e-18); ids=keep(1):keep(end); mass=mass(ids); origin=nodes(ids(1));
end

function [mass,origin]=convolveMass(a,x,b,y)
N=numel(a)+numel(b)-1; fftN=2^nextpow2(N);
mass=real(ifft(fft(a,fftN).*fft(b,fftN))); mass=mass(1:N); origin=x+y;
end
