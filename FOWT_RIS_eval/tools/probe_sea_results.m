root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'sea_wave'));
S=load(fullfile(root,'sea_wave','FigOutput_RIS_WindDirection','RIS_wind_direction_results.mat'));
n=numel(S.motions); U=zeros(n,1); dirs=zeros(n,1); tilt=zeros(n,1); heave=zeros(n,1);
for j=1:n
  q=S.motions{j};
  U(j)=q.parameters.windSpeed; dirs(j)=q.parameters.windAngleDeg;
  tilt(j)=hypot(std(q.rollDeg),std(q.pitchDeg));
  heave(j)=std(q.heave);
end
save(fullfile(root,'FOWT_RIS_eval','results','sea_model_envelope.mat'),'U','dirs','tilt','heave','-v7');
fprintf('saved %d cases\n',n);
