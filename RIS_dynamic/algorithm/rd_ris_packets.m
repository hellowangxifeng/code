function L = rd_ris_packets(c,g,packets)
%RD_RIS_PACKETS Lower-face MC packets -> LC-plane packets of the (moving) board.
%   packets is the second output of W_RIS/wris_mc for the posed geometry g.
%   Each packet is refracted water->glass (Snell + Fresnel), carried to the
%   LC screen at height g.hIn and kept only inside a 14.4 mm transparent tile
%   (same rounding and mask as RIS_phase/rp_illuminate). Packets keep their
%   raw MC power P0*w/N; no deterministic re-weighting is applied, so the
%   RIS and direct links share one unbiased estimator.
scale = c.source.power/packets.photons;
local = (packets.pos-g.pRis.')*g.RRis;
dirLocal = packets.dir*g.RRis;
qIn = c.optics.nWater*dirLocal(:,1:2);
kg = sqrt(max(c.optics.nGlass^2-sum(qIn.^2,2),realmin));
lc = local(:,1:2)+g.hIn*qIn./kg;
count = c.array.count;
pitch = c.array.pitchM;
ij = round(lc./pitch+(count+1)/2);
centres = (ij-(count+1)/2).*pitch;
use = all(ij>=1 & ij<=count,2) & all(abs(lc-centres)<=c.array.apertureM/2,2) & ...
    all(abs(lc)<=[c.ris.Lx c.ris.Ly]/2,2);
power = packets.weight(use)*scale.*rp_fresnel(c.optics.nWater,c.optics.nGlass,qIn(use,:));
L = struct('lcXY',lc(use,:),'qIn',qIn(use,:),'powerW',power, ...
    'unitId',ij(use,1)+(ij(use,2)-1)*count(1),'ballistic',packets.ballistic(use), ...
    'platePowerW',sum(packets.weight)*scale);
end
