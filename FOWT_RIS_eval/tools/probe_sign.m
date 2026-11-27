root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase'),fullfile(root,'FOWT_RIS_eval','src')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; q=s.cases{1}; g=q.geometry; a=q.illumination;
target.mode='direction'; target.u=[0 0 1];
p=fx_phase_design(c,g,a,target);
fprintf('cmd unique: '); disp(unique(p.commandDeg,'rows'));
d=mod(p.phaseAtReferenceRad-q.phase.phaseAtReferenceRad,2*pi); d=min(d,2*pi-d);
fprintf('piston diff vs cached: max %.3e rad, median %.3e\n',max(d),median(d));
% gradient diff
fprintf('grad diff max %.3e\n',max(abs(p.gradientOverK0(:)-q.phase.gradientOverK0(:))));
fprintf('outgoing valid diff %d\n',sum(p.outgoing.valid~=q.phase.outgoing.valid));
