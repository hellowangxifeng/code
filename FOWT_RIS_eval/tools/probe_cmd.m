root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; q=s.cases{1}; c=s.config;
cmd=q.phase.commandDeg;
lim=c.device.maxAngleDeg; step=c.device.angleStepDeg; maxc=floor(lim/step)*step;
fprintf('max command %.6f deg\n',maxc);
for j=1:2
  v=cmd(:,j);
  fprintf('axis %d: min %.4f max %.4f, at +lim %d, at -lim %d, of %d\n',j,min(v),max(v),sum(abs(v-maxc)<1e-9),sum(abs(v+maxc)<1e-9),numel(v));
  fprintf('  percentiles: '); fprintf('%.3f ',prctile(v,[0 5 25 50 75 95 100])); fprintf('\n');
end
% required command to hit +z for each unit
qin=q.phase.qInAtReference; 
fprintf('qin pct: \n'); disp(prctile(qin,[0 5 50 95 100]));
fprintf('b needed for +z = -qin; max |b| %.4f\n',c.optics.nAir*sind(maxc));
