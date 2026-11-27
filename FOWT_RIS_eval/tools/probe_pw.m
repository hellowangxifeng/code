root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; q=s.cases{1}; g=q.geometry; a=q.illumination;
ids=a.activeIds; ref=a.referenceXY(ids,:);
incoming=rp_layer_path(ref-g.sourceLocal(1:2).',[-g.sourceLocal(3),g.hIn],[c.optics.nWater,c.optics.nGlass],c.numerics.bisectionIterations);
qin=incoming.q; P=a.unitDirectPowerW(ids); P=P/sum(P); B=c.optics.nAir*sind(c.device.maxAngleDeg);
% candidate q grid in [-1,1]^2
qx=linspace(-1,1,801); qy=linspace(-1,1,801); [QX,QY]=meshgrid(qx,qy); RR=hypot(QX,QY);
% for each candidate, sum power of units with |q-qin|<=B axiswise
% do it by binning: use histcounts-like loop over coarse grid then refine
best=-1;
for ix=1:numel(qx)
  % vectorized over y: mask Nunits x Ny too big; loop y
end
% Instead: compute for each unit the feasible box center; power-weighted density on a coarse grid
nG=201;
gx=linspace(-1,1,nG); gy=gx;
dens=zeros(nG);
for i=1:numel(ids)
  lo=max(-1,qin(i,:)-B); hi=min(1,qin(i,:)+B);
  ix=discreteRange(lo(1),hi(1),gx); iy=discreteRange(lo(2),hi(2),gy);
  dens(iy,ix)=dens(iy,ix)+P(i);
end
[~,j]=max(dens(:)); [iy,ix]=ind2sub(size(dens),j);
fprintf('max power-covered q=[%.3f %.3f], angle %.2f deg, covered power %.3f\n',gx(ix),gy(iy),asind(hypot(gx(ix),gy(iy))),dens(j));
% weighted median-ish: choose q maximizing covered power along the mean direction
th=atan2d(mean(qin(:,2)),mean(qin(:,1)));
for ang=40:2:72
  qq=[cosd(th) sind(th)]*sind(ang);
  mask=all(abs(qin-qq)<=B,2) & sind(ang)<1;
  fprintf('angle %5.1f: covered power %.3f, units %d\n',ang,sum(P(mask)),sum(mask));
end

function idx=discreteRange(lo,hi,g)
  idx=find(g>=lo-1e-12 & g<=hi+1e-12);
end
