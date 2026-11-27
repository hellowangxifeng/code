root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase'),fullfile(root,'FOWT_RIS_eval','src')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; q=s.cases{1}; g=q.geometry; a=q.illumination; mc=q.mc;
I=load(fullfile(root,'FOWT_RIS_eval','results','improved_nominal.mat'));
target.mode='direction'; target.u=I.uWorld; target.curvature='cancel';
tic; p=fx_phase_design(c,g,a,target); fprintf('phase %.2f s\n',toc);
fprintf('angle mean %.3f std %.4f; controlled %.6f\n',mean(p.referenceAngleToWorldZDeg),std(p.referenceAngleToWorldZDeg),p.controlledPowerFraction);
r=rp_trace(c,g,a,p);
w=rp_receiver_window(c,g,a,p,r); p.centerXY=w.centerXY; p.receiverHalfWidthM=w.halfWidthM;
tic; wave=rp_wave(c,g,a,mc,p,r); fprintf('wave %.1f s\n',toc);
f=rp_footprint(c,g,p,wave);
[X,Y]=meshgrid(wave.x,wave.y); cx=f.centerXY(1); cy=f.centerXY(2);
fprintf('CURVATURE: center [%.4f %.4f] radius %.4f peak %.3f map %.5f\n',f.centerXY,f.radiusM,max(wave.irradianceWm2(:)),wave.mapPowerW);
for R=[0.01 0.025 0.05 0.1]
  m=(X-cx).^2+(Y-cy).^2<=R^2; P=sum(wave.irradianceWm2(m))*wave.cellAreaM2;
  fprintf('  R=%.3f P=%.5f (%.2f%%)\n',R,P,100*P/wave.mapPowerW);
end
save(fullfile(root,'FOWT_RIS_eval','results','improved_curvature.mat'),'p','wave','f','w','-v7.3');
