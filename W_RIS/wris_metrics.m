function s = wris_metrics(xy,w)
%WRIS_METRICS Weighted second moments; axes are full 4-sigma diameters.
% Halo-sensitive D4sigma axes describe the captured distribution, not a fit.
s=struct('center',[NaN NaN],'axes',[NaN NaN],'ratio',NaN,'orientationDeg',NaN);
if isempty(w) || sum(w)<=0, return; end
sw=sum(w); center=sum(xy.*w,1)/sw; q=xy-center;
covariance=q.'*(q.*w)/sw; [V,D]=eig(covariance,'vector');
[e,order]=sort(max(0,D),'descend'); diam=4*sqrt(e);
s.center=center; s.axes=diam.'; s.ratio=diam(1)/max(diam(2),realmin);
v=V(:,order(1)); s.orientationDeg=mod(atan2d(v(2),v(1))+90,180)-90;
end
