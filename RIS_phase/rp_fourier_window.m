function [F,fx,fy,df] = rp_fourier_window(values,spacing,padding,bounds,maximumPoints)
%RP_FOURIER_WINDOW Exact sampled Fourier integral on a selected frequency window.
% CZT evaluates the same uniformly spaced DFT as dense zero padding, without
% allocating the unneeded global spectrum. Samples are centred on the aperture.
arguments
    values (:,:) {mustBeNumeric,mustBeFinite}
    spacing (1,2) double {mustBeReal,mustBeFinite,mustBePositive}
    padding (1,2) double {mustBeReal,mustBeFinite,mustBeInteger,mustBePositive}
    bounds (2,2) double {mustBeReal,mustBeFinite}
    maximumPoints (1,1) double {mustBeInteger,mustBePositive} = 32000000
end
[ny,nx]=size(values); fullSize=padding.*[nx ny];
assert(all(mod(fullSize,2)==0),'rp:FourierSize','Equivalent DFT sizes must be even.');
assert(all(bounds(1,:)<=bounds(2,:)),'rp:FrequencyBounds','Frequency limits are reversed.');
df=1./(fullSize.*spacing);
low=max(-fullSize/2,floor(bounds(1,:)./df)-3);
high=min(fullSize/2-1,ceil(bounds(2,:)./df)+3);
if any(high<low)
    F=complex(zeros(0)); fx=[]; fy=[]; return
end
low=min(low,fullSize/2-4); high=max(high,low+3);
counts=high-low+1;
assert(prod(counts)<=maximumPoints,'rp:FourierWindow', ...
    'Receiver Fourier window needs %g points; increase wave.maximumFrequencySamples.',prod(counts));
fx=(low(1):high(1))*df(1); fy=(low(2):high(2))*df(2);
wx=exp(-2i*pi*df(1)*spacing(1)); ax=exp(2i*pi*fx(1)*spacing(1));
wy=exp(-2i*pi*df(2)*spacing(2)); ay=exp(2i*pi*fy(1)*spacing(2));
F=czt(values.',counts(1),wx,ax).';
F=czt(F,counts(2),wy,ay)*prod(spacing);
centre=([nx ny]-1).*spacing/2;
F=F.*exp(2i*pi*(fy(:)*centre(2)+fx*centre(1)));
end
