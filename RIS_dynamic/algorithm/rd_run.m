function R = rd_run(c)
%RD_RUN Frozen-sea RIS versus no-RIS study: compute, evaluate, save.
%   R = rd_run(c) with c from rd_config. For every wind U in
%   c.sea.windsMps one frozen instant is built (rd_scene); for every height
%   z in c.scene.distancesM all links are computed at the common receiver
%   (rd_condition). Each condition and each RIS design is a fingerprinted
%   checkpoint under c.cache.folder, so an interrupted run resumes and
%   post-processing changes never recompute optics. rd_evaluate selects the
%   schemes and builds every curve; the light result R is written to
%   <c.output>/ris_dynamic.mat together with CSV tables.
if nargin==0, c = rd_config(); end
rd_check(c);
if ~isfolder(c.output), mkdir(c.output); end
z = c.scene.distancesM;
U = c.sea.windsMps;
C = cell(numel(z),numel(U));
physical = rd_physical(c);
code = [rd_code_hash('direct') rd_code_hash('field')];
clock = tic;
for iw = 1:numel(U)
    files = arrayfun(@(zz) fullfile(c.cache.folder,sprintf('condition_U%g_z%g.mat',U(iw),zz)), ...
        z,'UniformOutput',false);
    prints = arrayfun(@(zz) struct('physical',physical,'wind',U(iw),'z',zz,'code',code), ...
        z,'UniformOutput',false);
    need = ~cellfun(@(f,p) rd_cache_valid(f,p,c.run.resume),files,prints);
    s = [];
    if any(need), s = rd_scene(c,U(iw)); end
    for iz = 1:numel(z)
        C{iz,iw} = rd_cached(files{iz},prints{iz},@() rd_condition(c,s,z(iz)),c.run.resume);
    end
    clear s
    fprintf('[run] U = %g m/s done (%.1f min)\n',U(iw),toc(clock)/60);
end
R = rd_evaluate(c,C);
save(fullfile(c.output,'ris_dynamic.mat'),'R','-v7.3');
rd_tables(R,c.output);
fprintf('\n[run] required P_t (dBm) for ABER = %g, main channel; * = best RIS optimizer\n',c.rank.targets(1));
for iw = 1:numel(U)
    for iz = 1:numel(z)
        k = R.cond(iz,iw).rank{1}.aber;
        mark = repmat({' '},1,5); mark(k.winners) = {'*'};
        fprintf('  U = %2g m/s, z = %2g m |',U(iw),z(iz));
        for j = 1:5
            fprintf(' %s %6.1f%s',R.schemeShort{j},k.requiredDBm(j),mark{j});
        end
        fprintf(' | clear: %d\n',k.clear);
    end
end
fprintf('[run] results: %s (%.1f min)\n',fullfile(c.output,'ris_dynamic.mat'),toc(clock)/60);
end
