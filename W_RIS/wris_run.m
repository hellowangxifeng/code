function s = wris_run(input)
%WRIS_RUN Calculate, validate and save; call wris_plot(s) separately to draw.
% s=wris_run(c) accepts the configuration edited in run_wris_irradiance.m.
% s=wris_run('smoke') uses 10000 photons per snapshot for a quick run.
if nargin==0, input='production'; end
wris_setup;
if isstruct(input), c=input; else, c=wris_config(input); end
if ~isfolder(c.output), mkdir(c.output); end
flotationAudit=audit_RIS_flotation(c.ris);
tic; s=wris_simulate(c); s.validation=wris_validate(c,s);
s.flotationAudit=flotationAudit; s.elapsedSeconds=toc;
save(fullfile(c.output,'wris_results.mat'),'s','-v7.3');
wris_report(s);
fprintf('Completed in %.1f seconds. Data: %s\n',s.elapsedSeconds,c.output);
end
