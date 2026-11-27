function R = rd_run(input)
%RD_RUN Compute and save the dynamic-sea RIS / rough-sea direct comparison.
%   R = rd_run()           production settings (rd_config)
%   R = rd_run('smoke')    tiny end-to-end check (results/smoke)
%   R = rd_run(c)          an edited rd_config structure
%
%   Pipeline (no figure is drawn here; call plot_ris_dynamic afterwards):
%     1. calm-sea RIS designs and receiver alignment per height (rd_design)
%     2. sea-aware width selection per wind on an independent sea (rd_train)
%     3. calm sea + every wind: paired frame-by-frame MC of all links
%        (rd_condition) -> instantaneous gains
%     4. incidence sweep, calm and one wind (optional)
%     5. ABER / OP / CDF / summaries (rd_evaluate)
%   Every stage is checkpointed under <output>/checkpoints and reused only
%   when its inputs and the physics code are unchanged.
%
%   Output (c.output): ris_dynamic.mat (variable R, read by the plots),
%   ris_dynamic_summary.csv, ris_dynamic_designs.csv, ris_dynamic_training.csv,
%   ris_dynamic_incidence.csv and run_log.txt.
if nargin==0, input = 'production'; end
rd_setup();
if isstruct(input)
    c = input;
else
    c = rd_config(input);
end
rd_check(c);
checkpoints = fullfile(c.output,'checkpoints');
folders = {c.output,c.cache.folder,checkpoints};
for j = 1:numel(folders)
    if ~isfolder(folders{j}), mkdir(folders{j}); end
end
diary(fullfile(c.output,'run_log.txt'));
closeLog = onCleanup(@() diary('off'));
fprintf('\n===== RIS_dynamic %s | profile %s | %s =====\n',c.version,c.profile, ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf('Output: %s\n',c.output);
timer = tic;
code = rd_code_hash();
phys = rd_physical(c);
file = @(name) fullfile(checkpoints,[name '.mat']);
resume = c.run.resume;
samplesStatic = rp_link_samples(c,c.study.samples,c.study.seed);
samplesDynamic = rp_link_samples(c,c.dynamic.samples,c.study.seed+1);
samplesTraining = rp_link_samples(c,c.training.samples,c.study.seed+2);
samplesIncidence = rp_link_samples(c,c.incidence.samples,c.study.seed+3);
z = c.scene.distancesM;
nz = numel(z);
incidence0 = c.scene.incidenceDeg;
winds = c.sea.windsMps;

%% 1. Calm-sea designs
designs = cell(1,nz);
for iz = 1:nz
    fp = struct('phys',phys,'code',code,'z',z(iz),'incidence',incidence0);
    designs{iz} = rd_cached(file(sprintf('design_z%04.1f_a%05.2f',z(iz),incidence0)),fp, ...
        @() rd_design(c,z(iz),incidence0,[],samplesStatic),resume);
end
designs = alignDirect(c,designs,max(winds),code,phys,file,'');

%% 2. Sea-aware selection
training = cell(1,numel(winds));
if c.training.enabled
    for iw = 1:numel(winds)
        fp = struct('phys',phys,'code',code,'wind',winds(iw),'designs',{designSignature(designs)});
        training{iw} = rd_cached(file(sprintf('training_U%g',winds(iw))),fp, ...
            @() rd_train(c,winds(iw),designs,samplesTraining),resume);
    end
end

%% 3. Calm sea and every wind
allWinds = [0 winds];
main = cell(1,numel(allWinds));
for iw = 1:numel(allWinds)
    seaAware = [];
    if allWinds(iw)>0 && c.training.enabled, seaAware = [training{iw-1}.chosen]; end
    S = rd_links(c,designs,seaAware,samplesDynamic);
    opts = struct('frameStepS',c.dynamic.frameStepS,'photons',c.dynamic.photons, ...
        'seedOffset',10000*iw,'calmFrames',c.dynamic.calmFrames,'maps',true, ...
        'snapshotTimeS',c.maps.snapshotTimeS,'label','main');
    fp = struct('phys',phys,'code',code,'wind',allWinds(iw),'links',S.signature,'opts',opts);
    fprintf('[condition] U = %g m/s\n',allWinds(iw));
    main{iw} = rd_cached(file(sprintf('main_U%g',allWinds(iw))),fp, ...
        @() rd_condition(c,allWinds(iw),S,opts),resume);
end

%% 4. Incidence sweep
incidence = struct('anglesDeg',[],'windMps',[],'blocks',{{}},'designs',{{}});
if c.incidence.enabled
    angles = c.incidence.anglesDeg;
    wind = c.incidence.windMps;
    iTrain = find(winds==wind,1);
    incidence.anglesDeg = angles;
    incidence.windMps = wind;
    incidence.blocks = cell(numel(angles),2);
    incidence.designs = cell(numel(angles),1);
    for ia = 1:numel(angles)
        da = cell(1,nz);
        for iz = 1:nz
            specs = nominalSpecs(designs{iz},training,iTrain,iz);
            fp = struct('phys',phys,'code',code,'z',z(iz),'incidence',angles(ia),'specs',specs);
            da{iz} = rd_cached(file(sprintf('incdesign_z%04.1f_a%05.2f',z(iz),angles(ia))),fp, ...
                @() rd_design(c,z(iz),angles(ia),specs,samplesIncidence),resume);
        end
        da = alignDirect(c,da,wind,code,phys,file,'inc');
        S = rd_links(c,da,[],samplesIncidence);
        conditionWinds = [0 wind];
        for k = 1:2
            opts = struct('frameStepS',c.incidence.frameStepS,'photons',c.incidence.photons, ...
                'seedOffset',500000+10000*ia+1000*k,'calmFrames',c.dynamic.calmFrames, ...
                'maps',false,'snapshotTimeS',c.maps.snapshotTimeS,'label', ...
                sprintf('incidence %g deg',angles(ia)));
            fp = struct('phys',phys,'code',code,'wind',conditionWinds(k),'links',S.signature, ...
                'opts',opts,'incidence',angles(ia));
            fprintf('[incidence] %g deg, U = %g m/s\n',angles(ia),conditionWinds(k));
            incidence.blocks{ia,k} = rd_cached(file(sprintf('incidence_a%05.2f_U%g', ...
                angles(ia),conditionWinds(k))),fp,@() rd_condition(c,conditionWinds(k),S,opts),resume);
        end
        incidence.designs{ia} = cellfun(@lightDesign,da,'UniformOutput',false);
    end
end

%% 5. Statistics and output
fprintf('[evaluate] ABER / OP / CDF\n');
R = rd_evaluate(c,designs,training,main,incidence);
R.codeHash = code;
R.elapsedSeconds = toc(timer);
R.createdAt = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
R.matlabVersion = version;
save(fullfile(c.output,'ris_dynamic.mat'),'R','-v7.3');
writetable(R.summary,fullfile(c.output,'ris_dynamic_summary.csv'));
writetable(R.designTable,fullfile(c.output,'ris_dynamic_designs.csv'));
writetable(R.trainingTable,fullfile(c.output,'ris_dynamic_training.csv'));
if ~isempty(R.incidenceTable)
    writetable(R.incidenceTable,fullfile(c.output,'ris_dynamic_incidence.csv'));
end
fprintf('Done in %.1f min. Results: %s\n',R.elapsedSeconds/60,fullfile(c.output,'ris_dynamic.mat'));
main1 = R.summary(R.summary.turbulenceId==1,:);
disp(main1(:,{'windMps','distanceM','scheme','meanGainDB','zeroFraction','aberFixed','opFixed'}));
end

function designs = alignDirect(c,designs,wind,code,phys,file,tag)
missing = find(cellfun(@(d) any(~isfinite(d.direct.center)),designs));
if isempty(missing) || ~(wind>0), return; end
g0 = designs{1}.g;
Z = cellfun(@(d) d.Z,designs);
fp = struct('phys',phys,'code',code,'wind',wind,'Z',Z,'pTx',g0.pTx,'dTx',g0.dTx);
centers = rd_cached(file(sprintf('%spilot_a%05.2f_U%g',tag,designs{1}.incidenceDeg,wind)),fp, ...
    @() rd_pilot(c,wind,g0,designs{1}.crossXY,Z),c.run.resume);
for iz = missing
    designs{iz}.direct.center = centers(iz,:);
    designs{iz}.direct.alignment = sprintf('rough-sea pilot, U = %g m/s',wind);
end
end

function specs = nominalSpecs(d,training,iTrain,iz)
% Keep the nominal family/width/margin of every RIS link (RIS_phase Fig10 rule).
keys = {d.ris.key};
specs = struct('key',keys,'family',{d.ris.family},'width',{d.ris.width}, ...
    'margin',{d.ris.margin});
seaAware = d.ris(strcmp(keys,'robust_aber'));
if ~isempty(iTrain) && ~isempty(training{iTrain}) && isfinite(training{iTrain}(iz).chosen)
    seaAware = d.bank(training{iTrain}(iz).chosen);
end
specs(end+1) = struct('key','sea_aware','family',seaAware.family,'width',seaAware.width, ...
    'margin',seaAware.margin);
end

function s = designSignature(designs)
s = cellfun(@(d) struct('family',{{d.bank.family}},'width',[d.bank.width], ...
    'margin',[d.bank.margin],'center',vertcat(d.bank.center)),designs,'UniformOutput',false);
end

function d = lightDesign(d)
d = rmfield(d,intersect(fieldnames(d),{'ris','bank'}));
end
