root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase'),fullfile(root,'FOWT_RIS_eval','src')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; q=s.cases{1}; g=q.geometry; a=q.illumination; mc=q.mc;
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
mask=dens>=0.999-1e-9; a2=ANG; a2(~mask)=Inf; [mn,j]=min(a2(:));
qstar=[QX(j),QY(j)];
fprintf('target q*=[%.4f %.4f], angle %.3f deg, coverage %.4f\n',qstar,mn,dens(j));
kz=sqrt(1-sum(qstar.^2)); uWorld=(g.RRis*[qstar,kz].').';
uWorld=uWorld/norm(uWorld);
fprintf('uWorld=[%.5f %.5f %.5f]\n',uWorld);
target.mode='direction'; target.u=uWorld;
tic; p=fx_phase_design(c,g,a,target); fprintf('phase design %.2f s\n',toc);
fprintf('angle mean %.3f, std %.4f, min %.3f max %.3f\n',mean(p.referenceAngleToWorldZDeg),std(p.referenceAngleToWorldZDeg),min(p.referenceAngleToWorldZDeg),max(p.referenceAngleToWorldZDeg));
fprintf('controlled power fraction %.6f\n',p.controlledPowerFraction);
fprintf('command at +lim x %d, y %d\n',sum(abs(p.commandDeg(:,1)+15.414)<1e-9),sum(abs(p.commandDeg(:,2)-15.414)<1e-9));
r=rp_trace(c,g,a,p);
w=rp_receiver_window(c,g,a,p,r); p.centerXY=w.centerXY; p.receiverHalfWidthM=w.halfWidthM;
tic; wave=rp_wave(c,g,a,mc,p,r); fprintf('wave %.1f s\n',toc);
f=rp_footprint(c,g,p,wave);
fprintf('IMPROVED: center [%.4f %.4f], radius %.4f m, threshold area %.5f, map power %.5f, peak %.3f\n', ...
  f.centerXY,f.radiusM,f.thresholdAreaM2,wave.mapPowerW,max(wave.irradianceWm2(:)));
fprintf('window halfwidth [%.4f %.4f]\n',p.receiverHalfWidthM);
save(fullfile(root,'FOWT_RIS_eval','results','improved_nominal.mat'),'p','wave','f','w','qstar','uWorld','-v7.3');
