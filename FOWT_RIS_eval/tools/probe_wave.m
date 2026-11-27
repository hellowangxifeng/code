root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s;
for k=1:numel(s.cases)
  q=s.cases{k};
  f=q.footprint; w=q.wave;
  fprintf(['case %d: angle mean %.3f | center [%.3f %.3f] centroid [%.3f %.3f] | ' ...
           'mapPower %.4f | rxP %.5g eval %d inWindow %d | pixel %.4f m\n'],k, ...
      mean(q.phase.referenceAngleToWorldZDeg),f.centerXY,f.centroidXY, ...
      w.mapPowerW,f.receivedPowerW,f.receivedPowerEvaluated,f.receiverApertureInWindow,f.pixelSizeM);
  fprintf('   map x [%.3f %.3f] y [%.3f %.3f] bins %dx%d, irradiance max %.4f mean %.5f\n', ...
      w.x(1),w.x(end),w.y(1),w.y(end),numel(w.x),numel(w.y),max(w.irradianceWm2(:)),mean(w.irradianceWm2(:)));
  fprintf('   radial r: '); fprintf('%.3f ',f.radialRadiusM); fprintf('\n   cover  : '); fprintf('%.3f ',f.radialCoverageFraction); fprintf('\n');
end
q=s.cases{1,1};
disp('RX position'); disp(s.rx{1}.position(:,1).');
disp('receiver window'); disp(q.receiverWindow);
disp('footprint rx'); fprintf('power %.5g evaluated %d radius %.4f\n',q.footprint.receivedPowerW,q.footprint.receivedPowerEvaluated,q.footprint.receiverApertureRadiusPixels*q.footprint.pixelSizeM);
