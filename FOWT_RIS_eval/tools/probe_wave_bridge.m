root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'FOWT_RIS_eval','src'));
G=load(fullfile(root,'FOWT_RIS_eval','data','processed','wave_gauge.mat'));
idx=8; t=G.(sprintf('t_%d',idx)); eta=G.(sprintf('eta_%d',idx));
w=fx_wave_from_gauge(t,eta(:,6),tStart=15,tEnd=35);
fprintf('%d components, Hs %.4f m, variance %.5f m^2\n',w.components,w.hsM,w.varianceM2);
% synthesise at x=0 over the segment and compare
sel=t>=15 & t<35; ts=t(sel)-t(find(sel,1)); e=eta(sel,6)-mean(eta(sel,6));
syn=real(exp(1i*(ts(:)*w.omega(:).'))*w.coefficient);
fprintf('synthesis rms error %.3e m vs std %.4f m\n',sqrt(mean((syn-e).^2)),std(e));
