function f = rd_spot(pk)
%RD_SPOT Power moments of the ballistic landing distribution (RIS_phase definitions).
%   equivalentRadiusM = sqrt(2 trace C) (= W for a uniform disc of radius W),
%   minor/majorRadiusM = 2 sqrt(eig C), r50M / r90M encircled-power radii.
%   If no ballistic packet arrives, all packets are used for the moments,
%   but ballisticPowerW stays 0 (used for receiver alignment decisions).
f = struct('powerW',0,'ballisticPowerW',0,'centerXY',[NaN NaN],'covarianceM2',nan(2), ...
    'equivalentRadiusM',NaN,'minorRadiusM',NaN,'majorRadiusM',NaN,'axisRatio',NaN, ...
    'r50M',NaN,'r90M',NaN);
if isempty(pk) || isempty(pk.power), return; end
f.powerW = sum(pk.power);
f.ballisticPowerW = sum(pk.power(pk.ballistic));
use = pk.ballistic;
if ~any(use), use = true(size(pk.power)); end
w = pk.power(use);
P = sum(w);
if ~(P>0), return; end
xy = pk.xy(use,:);
center = (w.'*xy)/P;
D = xy-center;
C = (D.*w).'*D/P;
C = (C+C.')/2;
axesM = 2*sqrt(max(eig(C),0));
[r,order] = sort(vecnorm(D,2,2));
F = cumsum(w(order))/P;
f.centerXY = center;
f.covarianceM2 = C;
f.equivalentRadiusM = sqrt(2*trace(C));
f.minorRadiusM = min(axesM);
f.majorRadiusM = max(axesM);
f.axisRatio = max(axesM)/max(min(axesM),realmin);
f.r50M = r(find(F>=0.5,1));
f.r90M = r(find(F>=0.9,1));
end
