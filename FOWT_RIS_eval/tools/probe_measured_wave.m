%FX_MEASURED_WAVE Probe: drive RIS_float with measured wave-gauge records.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'sea_wave'));
G=load(fullfile(root,'FOWT_RIS_eval','data','processed','wave_gauge.mat'));
for idx=[6 8]
  t=G.(sprintf('t_%d',idx)); eta=G.(sprintf('eta_%d',idx)); name=G.(sprintf('name_%d',idx));
  fs=1/median(diff(t)); sel=t>=15 & t<35;
  tt=t(sel)-t(find(sel,1)); e=eta(sel,6); e=e-mean(e);
  N=numel(e); T=N/fs; X=fft(e); f=(0:N-1).'/T; amp=2*abs(X)/N; ph=angle(X);
  keep=f>0.15 & f<3 & amp>2e-4;
  w=struct('omega',2*pi*f(keep),'kx',(2*pi*f(keep)).^2/9.81, ...
      'ky',zeros(nnz(keep),1),'coefficient',amp(keep).*exp(1i*ph(keep)));
  ts=(0:N-1).'/fs;
  syn=real(exp(1i*2*pi*(ts*f(keep).'))*w.coefficient);
  fprintf('record %s: fs %.2f, N %d, %d components, Hs %.3f m, synthesis rms error %.2e m\n', ...
      name,fs,N,nnz(keep),4*sqrt(sum(amp(keep).^2)/2),sqrt(mean((syn-e).^2)));
  p=get_default_RIS_params(); p.wave=w; p.duration=20; p.dt=0.01; p.startTime=0;
  p.initialization='waveSteady'; p.plot=false;
  tic; o=RIS_float(p); fprintf('  RIS_float %.1f s\n',toc);
  roll=o.rollDeg(:)-mean(o.rollDeg); pitch=o.pitchDeg(:)-mean(o.pitchDeg); heave=o.heave(:)-mean(o.heave);
  fprintf('  RIS model: heave std %.4f m, roll std %.3f deg, pitch std %.3f deg, max tilt %.3f deg\n', ...
      std(heave),std(roll),std(pitch),max(hypot(roll,pitch)));
  save(fullfile(root,'FOWT_RIS_eval','results',sprintf('measured_wave_ris_%d.mat',idx)),'w','o','e','tt','-v7');
end
