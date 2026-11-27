root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'sea_wave'));
p=get_default_RIS_params(); p.duration=1; p.dt=0.01; p.plot=false;
o=RIS_float(p);
f=o.hydroFrequencyHz(:); B=o.radiationDamping; A=o.radiationModel.infiniteAddedMass; M=o.massMatrix; K=o.stiffnessMatrix;
for j=1:3
  m=M(j,j)+A(j);
  T=2*pi/sqrt(K(j,j)/m); fn=1/T;
  bj=interp1(f,B(j,:),fn,'linear','extrap');
  zeta=bj/(2*sqrt(m*K(j,j)));
  fprintf('mode %d: Tn %.3f s fn %.3f Hz, B %.4g, zeta %.4f\n',j,T,fn,bj,zeta);
end
fprintf('B heave at 0.4/0.8/1.2 Hz: '); fprintf('%.4g ',interp1(f,B(1,:),[0.4 0.8 1.2])); fprintf('\n');
