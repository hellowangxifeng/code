root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase')); rp_setup();
M=load(fullfile(root,'FOWT_RIS_eval','results','current_quadratic.mat'));
w=M.wave; f=M.f; [X,Y]=meshgrid(w.x,w.y);
fprintf('center %.4f %.4f, radius %.4f, area %.5f, cov %.3f\n',f.centerXY,f.radiusM,f.thresholdAreaM2,f.coverageFraction);
fprintf('radial r: '); fprintf('%.3f ',f.radialRadiusM(1:10:end)); fprintf('\n');
fprintf('cover  : '); fprintf('%.3f ',f.radialCoverageFraction(1:10:end)); fprintf('\n');
[~,im]=max(w.irradianceWm2(:)); [iy,ix]=ind2sub(size(w.irradianceWm2),im);
fprintf('peak at x %.4f y %.4f val %.3f\n',w.x(ix),w.y(iy),w.irradianceWm2(im));
for R=[0.005 0.01 0.02 0.03 0.05]
  m=(X-f.centerXY(1)).^2+(Y-f.centerXY(2)).^2<=R^2;
  fprintf('R %.3f: sum %.5f\n',R,sum(w.irradianceWm2(m))*w.cellAreaM2);
end
% sum of map
fprintf('map %.5f; total pixels %.2f\n',w.mapPowerW,numel(w.x)*numel(w.y));
