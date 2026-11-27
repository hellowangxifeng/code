function [local,range,inside] = wris_intersect(p,d,center,R,boardSize)
%WRIS_INTERSECT World rays -> body plane; row vectors use (xW-pW)*R.
% Active body-to-world R=Rz(yaw)Ry(pitch)Rx(roll), reused RIS convention.
normal = R(:,3);
den = d*normal;
range = ((center.'-p)*normal)./den;
valid = abs(den)>1e-12 & range>=0;
range(~valid) = Inf;
local = (p+d.*range-center.')*R;
local(~valid,:) = NaN;
inside = valid & abs(local(:,1))<=boardSize(1)/2 & abs(local(:,2))<=boardSize(2)/2;
end
