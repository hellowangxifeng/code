root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'sea_wave'));
S=load(fullfile(root,'sea_wave','FigOutput_RIS_WindDirection','RIS_wind_direction_results.mat'));
m=S.motions{1};
f=m.hydroFrequencyHz(:); R=m.responseAmplitudeOperator6;
disp(class(f)); disp(size(f)); disp(class(R)); disp(size(R));
disp(f(1:3)); disp(R(3,1:3));
