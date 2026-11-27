root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'sea_wave'));
S=load(fullfile(root,'sea_wave','FigOutput_RIS_WindDirection','RIS_wind_direction_results.mat'));
m=S.motions{1}; f=m.hydroFrequencyHz(:); R=m.responseAmplitudeOperator6;
names={'surge','sway','heave','roll','pitch','yaw'};
for j=1:6
  [~,i]=max(R(j,:));
  fprintf('%-6s peak %.4f Hz  value %.4g\n',names{j},f(i),R(j,i));
end
for ff=[0.4 0.8 1.2]
  [~,i]=min(abs(f-ff));
  fprintf('f=%.2f Hz: heave %.3f m/m, pitch %.3f deg/m\n',ff,R(3,i),R(5,i));
end
save(fullfile(root,'FOWT_RIS_eval','results','sea_model_rao.mat'),'f','R','names','-v7');
