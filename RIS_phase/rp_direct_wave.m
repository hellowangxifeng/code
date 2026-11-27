function b=rp_direct_wave(c,g)
%RP_DIRECT_WAVE Continuous water-air direct beam without glass or LC masks.
% Rotated horizontal grids resolve the extreme astigmatism near critical
% incidence; the receiver calculation transforms the SAME physical samples.
folder=fullfile(c.output,'wave_cache');
file=fullfile(folder,sprintf('direct_a%05.2f_z%04.1f.mat',g.incidenceDeg, ...
    g.rxPlaneZ-g.pRis(3)-c.ris.thickness));
modelSource=[fileread([mfilename('fullpath') '.m']),fileread(which('rp_direct_density')), ...
    fileread(which('rp_direct_channel'))];
if isfile(file)
    old=load(file,'b','signature','model');
    if strcmp(old.signature,c.study.signature) && strcmp(old.model,modelSource)
        b=old.b; return
    end
end
az=atan2(g.dTx(2),g.dTx(1)); Q=[cos(az),-sin(az);sin(az),cos(az)];
slope=linspace(-.08,.08,641); ds=slope(2)-slope(1);
[u,v]=meshgrid(slope); directions=[u(:),v(:),ones(numel(u),1)];
directions=directions./vecnorm(directions,2,2)*g.RTx.';
q=c.optics.nWater*directions(:,1:2); kw=c.optics.nWater*directions(:,3);
ka=sqrt(max(0,1-sum(q.^2,2))); depth=-g.sourceLocal(3);
height=g.rxPlaneZ-g.pRis(3); good=kw>0 & ka>0;
power=2*c.source.power/(pi*c.source.divergence^2)* ...
    exp(-2*(u(:).^2+v(:).^2)/c.source.divergence^2)*ds^2;
power=power.*exp(-(c.water.absorption+c.water.scattering)*depth./max(directions(:,3),realmin));
power=power.*rp_fresnel(c.optics.nWater,c.optics.nAir,q); power(~good)=0; total=sum(power);
if total<1e-100
    center=g.pRis(1:2).'*Q; half=[.4 .4]; count=[192 192];
else
    worldXY=g.pRis(1:2).'+g.sourceLocal(1:2).'+depth*q(good,:)./kw(good)+height*q(good,:)./ka(good);
    xy=worldXY*Q; weight=power(good); lo=zeros(1,2); hi=lo; core=lo;
    for dimension=1:2
        lo(dimension)=quantileWeighted(xy(:,dimension),weight,5e-5);
        hi(dimension)=quantileWeighted(xy(:,dimension),weight,1-5e-5);
        core(dimension)=quantileWeighted(xy(:,dimension),weight,.9)-quantileWeighted(xy(:,dimension),weight,.1);
    end
    center=(lo+hi)/2; half=max((hi-lo)/2*1.05,[.4 .4]);
    count=max([192 192],min([4096 512],ceil(160*half./max(core,1e-6))));
end
x=linspace(center(1)-half(1),center(1)+half(1),count(1));
y=linspace(center(2)-half(2),center(2)+half(2),count(2));
if total<1e-100
    I=zeros(count(2),count(1));
else
    [X,Y]=meshgrid(x,y); points=[X(:),Y(:)]*Q.';
    I=reshape(rp_direct_density(c,g,points),count(2),count(1));
end
w=struct('x',x,'y',y,'irradianceWm2',I,'cellAreaM2',(x(2)-x(1))*(y(2)-y(1)), ...
    'model','continuous Gaussian angular density, exact Snell ray Jacobian; direct beam only', ...
    'worldToMap',Q,'angularCubaturePowerW',total,'scatteredIrradianceWm2',zeros(size(I)));
f=rp_spot_stats(w); if total<1e-100, f.centerXY=center; end
f.centerXY=f.centerXY*Q.'; f.covarianceM2=Q*f.covarianceM2*Q.';
h=rp_direct_channel(c,g,w,f.centerXY,rp_link_samples(c));
b=struct('family','direct','targetWidthM',0,'wave',w,'spot',f,'channel',h, ...
    'designMetrics',rp_link_metrics(c,h,c.study.designPowerW));
model=modelSource; signature=c.study.signature; save(file,'b','signature','model','-v7.3');
end

function q=quantileWeighted(x,w,p)
[x,order]=sort(x); F=cumsum(w(order)); q=x(find(F>=p*F(end),1));
end
