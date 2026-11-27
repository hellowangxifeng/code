function w=rp_proxy_wave(c,g,a,p)
%RP_PROXY_WAVE Screening only: deterministic aperture rays plus diffraction blur.
N=c.study.raySamplesPerAxis;
axis=((1:N)-0.5)/N-0.5;
[X,Y]=meshgrid(axis*c.array.apertureM(1),axis*c.array.apertureM(2));
n=numel(p.unitIds); count=numel(X);
xy=repelem(a.unitCenters(p.unitIds,:),count,1)+repmat([X(:),Y(:)],n,1);
input=rp_incident_field(c,g,xy);
packets=struct('lcXY',xy,'qIn',input.q, ...
    'sampleUnitIndex',repelem((1:n).',count),'unitId',repelem(p.unitIds,count), ...
    'powerW',input.irradianceWm2*prod(c.array.apertureM)/count);
r=rp_trace(c,g,packets,p);
use=r.valid & r.powerW>max(r.powerW)*1e-10;
points=r.landing(use,1:2); power=r.powerW(use);
center=sum(points.*power,1)/sum(power);
half=max(max(abs(points-center),[],1)+0.015,[0.30 0.30]);
bins=c.study.proxyBins;
xe=linspace(center(1)-half(1),center(1)+half(1),bins+1);
ye=linspace(center(2)-half(2),center(2)+half(2),bins+1);
ix=discretize(points(:,1),xe); iy=discretize(points(:,2),ye);
area=(xe(2)-xe(1))*(ye(2)-ye(1));
I=accumarray([iy ix],power,[bins bins],@sum,0)/area;
% A width of at least one pixel prevents ray delta spikes from winning screening.
sigma=max(1,2*c.source.wavelength*(g.rxPlaneZ-g.pRis(3))/min(c.array.apertureM)/sqrt(area));
t=-ceil(4*sigma):ceil(4*sigma); h=exp(-t.^2/(2*sigma^2)); h=h/sum(h);
I=conv2(h,h,I,'same');
w=struct('x',(xe(1:end-1)+xe(2:end))/2,'y',(ye(1:end-1)+ye(2:end))/2, ...
    'irradianceWm2',I,'cellAreaM2',area,'model','screening ray envelope; NOT coherent wave');
end
