root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; q=s.cases{1}; g=q.geometry; a=q.illumination;
ids=a.activeIds; ref=a.referenceXY(ids,:);
incoming=rp_layer_path(ref-g.sourceLocal(1:2).',[-g.sourceLocal(3),g.hIn],[c.optics.nWater,c.optics.nGlass],c.numerics.bisectionIterations);
qin=incoming.q; P=a.unitDirectPowerW(ids); P=P/sum(P); B=c.optics.nAir*sind(c.device.maxAngleDeg);
nG=401; gx=linspace(-1,1,nG); gy=gx; dens=zeros(nG);
for i=1:numel(ids)
  lo=max(-1,qin(i,:)-B); hi=min(1,qin(i,:)+B);
  ix=find(gx>=lo(1)-1e-12 & gx<=hi(1)+1e-12); iy=find(gy>=lo(2)-1e-12 & gy<=hi(2)+1e-12);
  dens(iy,ix)=dens(iy,ix)+P(i);
end
[QX,QY]=meshgrid(gx,gy); ANG=asind(min(1,hypot(QX,QY)));
for target=[0.95 0.99 0.995 0.999 1.0]
  mask=dens>=target-1e-9;
  if ~any(mask(:)), fprintf('target %.3f infeasible\n',target); continue; end
  a2=ANG; a2(~mask)=Inf; [mn,j]=min(a2(:));
  fprintf('coverage>=%.3f: min angle %.2f deg at q=[%.3f %.3f]; max coverage at <=50/55/60 deg: ',target,mn,QX(j),QY(j));
  for alim=[50 55 60]
    m2=mask & ANG<=alim; fprintf('%.4f ',max(dens(m2)));
  end
  fprintf('\n');
end
% Pareto: min angle for each coverage level
fprintf('\ncoverage at fixed angles (best azimuth):\n');
for alim=45:1:65
  m2=ANG<=alim; fprintf('%d: %.4f  ',alim,max(dens(m2)));
  if mod(alim,5)==0, fprintf('\n'); end
end
fprintf('\n');
