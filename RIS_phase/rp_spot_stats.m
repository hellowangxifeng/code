function f=rp_spot_stats(w)
%RP_SPOT_STATS Power moments and energy radii, without a peak threshold.
[X,Y]=meshgrid(w.x,w.y); I=w.irradianceWm2;
P=sum(I,'all')*w.cellAreaM2;
weight=I(:)/max(sum(I,'all'),realmin);
center=[sum(X(:).*weight),sum(Y(:).*weight)];
D=[X(:),Y(:)]-center; C=(D.*weight).'*D;
axesM=2*sqrt(max(eig((C+C.')/2),0));
[radius,idx]=sort(vecnorm(D,2,2)); F=cumsum(weight(idx));
r50=radius(find(F>=0.5,1)); r90=radius(find(F>=0.9,1));
if isempty(r50), r50=NaN; r90=NaN; end
f=struct('powerW',P,'centerXY',center,'covarianceM2',C, ...
    'equivalentRadiusM',sqrt(2*trace(C)),'minorRadiusM',min(axesM), ...
    'majorRadiusM',max(axesM),'axisRatio',max(axesM)/max(min(axesM),realmin), ...
    'r50M',r50,'r90M',r90,'concentrationR50R90',r50/max(r90,realmin), ...
    'peakWm2',max(I,[],'all'));
end
