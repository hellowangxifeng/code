root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'FOWT_RIS_eval','src'),fullfile(root,'sea_wave'),fullfile(root,'RIS_phase'));
raw=fullfile(root,'浮式平台六自由度实验数据','FOWT Platform Testing','FOWTmotion_11_3_45mmp8Hz.mat');
m=fx_motion_from_experiment(raw);
fprintf('t %d samples, %.2f s, fs %.1f; bottom z range %.4f..%.4f m\n',numel(m.t),m.t(end),1/median(diff(m.t)),min(m.bottomCenter(3,:)),max(m.bottomCenter(3,:)));
fprintf('mean roll/pitch/yaw %.3f/%.3f/%.3f deg\n',mean(m.rollDeg),mean(m.pitchDeg),mean(m.yawDeg));
for k=[1 100 500]
  R=m.rotationMatrix(:,:,k); fprintf('orthogonality error %.2e, det %.6f\n',norm(R.'*R-eye(3)),det(R));
end
fprintf('normal z range %.4f..%.4f\n',min(m.normal(3,:)),max(m.normal(3,:)));
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
c=S.s.config;
tx=struct('position',repmat([-5;7.0711;-5.0034],1,numel(m.t)), ...
    'direction',repmat([0.5;-0.7071;0.5]/norm([0.5;-0.7071;0.5]),1,numel(m.t)), ...
    'rotation',repmat(eye(3),1,1,numel(m.t)));
rx=struct('position',repmat([0;0;5],1,numel(m.t)),'nominal',[0;0;5]);
g=rp_geometry(c,m,tx,rx,1000);
fprintf('rp_geometry works: pRis [%.4f %.4f %.4f], incidence %.3f deg\n',g.pRis,g.incidenceDeg);
