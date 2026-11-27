function map=rp_fig1_mc_map(packets)
%RP_FIG1_MC_MAP MC density estimate and full-packet wide-field context map.
arguments
    packets (1,1) struct
end
p=packets.powerW; xy=packets.landingXY; core=packets.ballistic;
assert(any(core),'rp:DirectCore','Fig1 requires a transmitted ballistic main beam.');
center=sum(xy(core,:).*p(core),1)/sum(p(core));
D=xy(core,:)-center; C=(D.*p(core)).'*D/sum(p(core));
effectiveN=sum(p(core))^2/sum(p(core).^2);
% Scott bandwidth for a two-dimensional weighted MC density estimator.
bandwidth=sqrt(diag(C)).'*effectiveN^(-1/6);
halfWidth=4*sqrt(diag(C)).'+6*bandwidth;
count=max(128,ceil(2*halfWidth./(bandwidth/2))); count=2*ceil(count/2);
edgesX=linspace(center(1)-halfWidth(1),center(1)+halfWidth(1),count(1)+1);
edgesY=linspace(center(2)-halfWidth(2),center(2)+halfWidth(2),count(2)+1);
raw=histogramPower(xy,p,edgesX,edgesY); dx=diff(edgesX(1:2));dy=diff(edgesY(1:2));
I=smoothDensity(raw/(dx*dy),bandwidth./[dx dy]);
Ihalf=smoothDensity(raw/(dx*dy),bandwidth./[dx dy]/2);
w=struct('x',(edgesX(1:end-1)+edgesX(2:end))/2, ...
    'y',(edgesY(1:end-1)+edgesY(2:end))/2,'irradianceWm2',I,'cellAreaM2',dx*dy);
% The main panel resolves the beam; invertible asinh coordinates in the inset
% show EVERY transmitted MC packet, including grazing-angle diffuse tails.
offset=xy-center; scale=1; transformed=asinh(offset/scale);
extent=max(abs(transformed),[],1)*1.03+1e-6;
u=linspace(-extent(1),extent(1),161);v=linspace(-extent(2),extent(2),161);
haloPower=histogramPower(transformed,p,u,v);
physicalArea=diff(scale*sinh(v)).'*diff(scale*sinh(u));
map=struct('wave',w,'spot',rp_spot_stats(w),'rawPixelPowerW',raw, ...
    'bandwidthM',bandwidth,'halfBandwidthPeakWm2',max(Ihalf,[],'all'), ...
    'mainBoundsXY',[edgesX([1 end]);edgesY([1 end])], ...
    'mainPacketPowerW',sum(raw,'all'),'mainMapPowerW',sum(I,'all')*dx*dy, ...
    'ballisticPowerFractionInMainView',sum(p(core & all(abs(offset)<=halfWidth,2)))/sum(p(core)), ...
    'totalPowerFractionInMainView',sum(raw,'all')/sum(p), ...
    'halo',struct('u',(u(1:end-1)+u(2:end))/2,'v',(v(1:end-1)+v(2:end))/2, ...
        'irradianceWm2',haloPower./physicalArea,'pixelPowerW',haloPower, ...
        'scaleM',scale,'centerXY',center,'capturedPowerFraction',sum(haloPower,'all')/sum(p)));
end

function mass=histogramPower(xy,p,xe,ye)
ix=discretize(xy(:,1),xe);iy=discretize(xy(:,2),ye);valid=isfinite(ix)&isfinite(iy);
mass=accumarray([iy(valid),ix(valid)],p(valid),[numel(ye)-1 numel(xe)-1],@sum,0);
end

function density=smoothDensity(raw,sigmaPixels)
ax=-ceil(5*sigmaPixels(1)):ceil(5*sigmaPixels(1));
ay=-ceil(5*sigmaPixels(2)):ceil(5*sigmaPixels(2));
kx=exp(-.5*(ax/sigmaPixels(1)).^2);ky=exp(-.5*(ay/sigmaPixels(2)).^2);
density=conv2(ky.'/sum(ky),kx/sum(kx),raw,'same');
end
