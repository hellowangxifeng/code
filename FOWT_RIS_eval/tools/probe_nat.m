root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'sea_wave'));
p=get_default_RIS_params(); p.duration=1; p.dt=0.01; p.plot=false;
o=RIS_float(p);
disp('massMatrix'); disp(o.massMatrix);
disp('stiffnessMatrix'); disp(o.stiffnessMatrix);
M=o.massMatrix; K=o.stiffnessMatrix;
for j=1:3
  fprintf('mode %d: Tn %.4f s (f %.4f Hz)\n',j,2*pi/sqrt(K(j,j)/M(j,j)),sqrt(K(j,j)/M(j,j))/(2*pi));
end
fprintf('Ainf %s\n',mat2str(o.radiationModel.infiniteAddedMass(:).',4));
