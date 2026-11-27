function m = rd_link_metrics(c,G,psi,powerW,ocean,atmosphere)
%RD_LINK_METRICS ABER, OP and mean power with the exact FOV tail (RIS_phase eqs. 99-100).
%   G(s)   instantaneous aperture gain of receiver sample s (no fading),
%   psi(s) pointing error of the same sample. With the exact Rice tail
%   p_F = P(psi > Theta) and the in-FOV samples A = {psi <= Theta}:
%     ABER = p_F/2 + (1-p_F) E[ E_H[ Q(R P_t G H / 2 sigma_n) ] | A ]
%     OP   = p_F   + (1-p_F) E[ F_H( P_th / (P_t G) ) | A ]
%     <P_r> = (1-p_F) P_t E[G | A]
%   H = H_o H_a is the mean-one ocean EGG x atmospheric Gamma-Gamma fading
%   named by (ocean, atmosphere); the inner expectations are
%   RIS_phase/rp_link_metrics (log-domain convolution, thermal-noise OOK).
cc = c;
cc.study.turbulence = ocean;
cc.study.atmosphere = atmosphere;
theta = c.study.fovHalfAngleRad;
if c.study.angleSigmaRad>0
    pF = ncx2cdf((theta/c.study.angleSigmaRad)^2,2, ...
        sum(c.study.angleBiasRad.^2)/c.study.angleSigmaRad^2,'upper');
else
    pF = double(norm(c.study.angleBiasRad)>theta);
end
accepted = double(psi(:))<=theta;
assert(any(accepted),'rd:FOV','No receiver sample lies inside the FOV.');
h = struct('gain',double(G(accepted)),'fovOutageExact',0,'fovOutageSample',0);
r = rp_link_metrics(cc,h,powerW,c.metrics.order);
m = struct('aber',(1-pF)*r.aber+pF/2,'op',(1-pF)*r.op+pF, ...
    'meanPowerW',(1-pF)*r.meanReceivedPowerW,'meanGain',(1-pF)*mean(h.gain), ...
    'conditionalGain',mean(h.gain),'fovOutageExact',pF);
end
