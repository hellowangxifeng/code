function out = wrd_score(c,hits)
%WRD_SCORE Estimate aperture power, uncertainty and common-plane irradiance.
weights = hits(:,3);
aperture = sum(hits(:,1:2).^2,2)<=c.receiver.radius^2;
accepted = aperture & hits(:,5)>=cosd(c.receiver.fovDeg);
ballistic = hits(:,4)>0;
n = c.mc.photons;
scale = c.source.power/n;
powerSum = sum(weights(accepted));
out.rxPowerW = powerSum*scale;
out.ballisticPowerW = sum(weights(accepted & ballistic))*scale;
out.scatterPowerW = sum(weights(accepted & ~ballistic))*scale;
out.rxOnlyApertureW = sum(weights(aperture))*scale;
out.hitCount = nnz(accepted);
out.powerSEW = c.source.power*sqrt(max(0,sum(weights(accepted).^2)-powerSum^2/n)/(n*(n-1)));
edges = c.receiver.mapEdges;
out.axisM = (edges(1:end-1)+edges(2:end))/2;
ix = discretize(hits(:,1),edges);
iy = discretize(hits(:,2),edges);
inMap = isfinite(ix) & isfinite(iy);
shape = repmat(numel(edges)-1,1,2);
out.mapWm2 = accumarray([iy(inMap),ix(inMap)],weights(inMap),shape)*scale/diff(edges(1:2))^2;
out.mapPowerW = sum(weights(inMap))*scale;
end
