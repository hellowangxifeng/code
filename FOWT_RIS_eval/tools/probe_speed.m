root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; q=s.cases{1}; g=q.geometry; a=q.illumination; mc=q.mc; p=q.phase; r=q.rays;
tic; w=rp_wave(c,g,a,mc,p,r); fprintf('rp_wave elapsed %.2f s\n',toc);
fprintf('map power %.6f vs %.6f\n',w.mapPowerW,q.wave.mapPowerW);
