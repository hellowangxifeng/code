root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config;
I=load(fullfile(root,'FOWT_RIS_eval','results','improved_nominal.mat'));
for which=1:2
  if which==1, w=s.cases{1}.wave; f=s.cases{1}.footprint; tag='current'; else, w=I.wave; f=I.f; tag='improved'; end
  [XX,YY]=meshgrid(w.x,w.y);
  cx=f.centerXY(1); cy=f.centerXY(2);
  for R=[0.01 0.025 0.05 0.1]
    m=(XX-cx).^2+(YY-cy).^2<=R^2;
    P=sum(w.irradianceWm2(m))*w.cellAreaM2;
    fprintf('%s aperture R=%.3f: power %.5f W (%.2f%% of map)\n',tag,R,P,100*P/w.mapPowerW);
  end
  rr=hypot(XX-cx,YY-cy); [rs,idx]=sort(rr(:)); cw=cumsum(w.irradianceWm2(idx))*w.cellAreaM2; cw=cw/cw(end);
  [ru,iu]=unique(rs); cwu=cw(iu);
  for R=[0.01 0.025 0.05 0.1 0.2]
    fprintf('   %s radial power within %.3f m: %.4f\n',tag,R,interp1(ru,cwu,R));
  end
end
[~,im]=max(I.wave.irradianceWm2(:)); [iy,ix]=ind2sub(size(I.wave.irradianceWm2),im);
fprintf('improved peak at x=%.4f y=%.4f, center %.4f %.4f\n',I.wave.x(ix),I.wave.y(iy),I.f.centerXY);
