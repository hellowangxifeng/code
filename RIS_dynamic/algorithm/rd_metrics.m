function m = rd_metrics(c,G,rx,powerW)
%RD_METRICS Frame-wise and time-averaged mean power, ABER and OP.
%   G is samples x frames of instantaneous geometric gains (no fading).
%   Each frame is passed to RIS_phase/rp_link_metrics (joint mean-one EGG x
%   Gamma-Gamma fading, thermal-noise OOK, SNR outage) and the frames are
%   averaged with equal weight, i.e. time is one more uniformly sampled
%   random variable: ABER = E_t E_jitter E_h[BER], OP = E_t P(SNR < th | t).
G = double(G);
powerW = powerW(:).';
nf = size(G,2);
np = numel(powerW);
m = struct('powerW',powerW,'aberFrames',nan(nf,np),'opFrames',nan(nf,np), ...
    'meanGainFrames',nan(nf,1),'aber',nan(1,np),'op',nan(1,np),'meanGain',NaN);
if isempty(G) || any(~isfinite(G(:))), return; end
for f = 1:nf
    h = struct('gain',G(:,f),'fovOutageExact',rx.fovOutageExact, ...
        'fovOutageSample',rx.fovOutageSample);
    r = rp_link_metrics(c,h,powerW);
    m.aberFrames(f,:) = r.aber;
    m.opFrames(f,:) = r.op;
end
m.meanGainFrames = mean(G,1).';
m.aber = mean(m.aberFrames,1);
m.op = mean(m.opFrames,1);
m.meanGain = mean(m.meanGainFrames);
end
