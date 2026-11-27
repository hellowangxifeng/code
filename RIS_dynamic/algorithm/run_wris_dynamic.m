function results = run_wris_dynamic(input)
%RUN_WRIS_DYNAMIC Compute and save the paired rough-sea direct / moving-RIS links.
%   results = run_wris_dynamic()          production settings
%   results = run_wris_dynamic('smoke')   small end-to-end check
%   results = run_wris_dynamic(c)         an edited WRD_CONFIG structure
%
%   For every wind condition one directional sea is generated. The SAME wave
%   components drive the RIS six-DOF motion (sea_wave/RIS_float) and the rough
%   water-air boundary crossed by the direct link. At every frame all links
%   share the TX pose, beam, water, receiver and MC seed:
%     direct : TX -> water (HG MC) -> rough sea surface (Snell/Fresnel/TIR) -> air -> RX
%     flat   : the same, with the still-water plane z = 0 (old W_RIS baseline)
%     ris    : TX -> water (W_RIS MC) -> moving RIS (LC steering) -> air -> RX
%   No figure is drawn here; call plot_wris_dynamic afterwards.
%
%   Output (c.output): wris_dynamic.mat, wris_dynamic_frames.csv,
%   wris_dynamic_cases.csv, run_log.txt and checkpoints/condition_*.mat.
if nargin==0, input = 'production'; end
wrd_setup();
if isstruct(input)
    c = input;
else
    c = wrd_config(input);
end
wrd_check_dependencies();
checkConfig(c);
if ~isfolder(c.output), mkdir(c.output); end
checkpoints = c.run.checkpointFolder;
if isempty(checkpoints), checkpoints = fullfile(c.output,'checkpoints'); end
if ~isfolder(checkpoints), mkdir(checkpoints); end
logFile = fullfile(c.output,'run_log.txt');
diary(logFile);
closeLog = onCleanup(@() diary('off'));
fprintf('\n===== RIS_dynamic %s | profile %s | %s =====\n',c.version,c.profile, ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf('Output: %s\n',c.output);

conditions = c.comparison.conditions;
nRows = size(c.comparison.rotationXYDeg,1);
cases = cell(nRows,size(conditions,1));
codeHash = algorithmHash();
runTimer = tic;
for seaId = 1:size(conditions,1)
    file = fullfile(checkpoints,sprintf('condition_%02d_U%g_dir%g.mat',seaId, ...
        conditions(seaId,1),conditions(seaId,2)));
    fingerprint = struct('config',physicalConfig(c),'condition',conditions(seaId,:), ...
        'codeHash',codeHash);
    block = [];
    if c.run.resume && isfile(file)
        saved = load(file,'block','fingerprint');
        if isequaln(saved.fingerprint,fingerprint)
            block = saved.block;
            fprintf('Wind %g m/s: reusing checkpoint %s\n',conditions(seaId,1),file);
        else
            fprintf('Wind %g m/s: checkpoint does not match the current config/code, recomputing.\n', ...
                conditions(seaId,1));
        end
    end
    if isempty(block)
        block = runCondition(c,seaId,runTimer);
        save(file,'block','fingerprint','-v7');
    end
    cases(:,seaId) = block;
end

[frameTable,caseTable] = wrd_summary(cases,c);
results = struct('config',c,'cases',{cases},'frameTable',frameTable, ...
    'caseTable',caseTable,'codeHash',codeHash,'elapsedSeconds',toc(runTimer), ...
    'createdAt',char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')), ...
    'matlabVersion',version);
save(fullfile(c.output,'wris_dynamic.mat'),'results','-v7');
writetable(frameTable,fullfile(c.output,'wris_dynamic_frames.csv'));
writetable(caseTable,fullfile(c.output,'wris_dynamic_cases.csv'));
fprintf('Saved %d cases (%d frames each) in %.1f min to %s\n',numel(cases), ...
    numel(cases{1}.frames),results.elapsedSeconds/60,c.output);
disp(caseTable(:,{'windMps','nominalIncidenceDeg','link','meanW','outageProbability','zeroHitFraction'}));
end

function block = runCondition(c,seaId,runTimer)
cfg = c;
cfg.ris.windSpeed = c.comparison.conditions(seaId,1);
cfg.ris.windAngleDeg = c.comparison.conditions(seaId,2);
seaDefault = get_default_sea_params();
wave = RIS_wave_components(cfg.ris,cfg.sea.quadrature);
wave.coefficient = wave.coefficient*sqrt(cfg.sea.neumannConstant/seaDefault.neumannConstant);
cfg.ris.wave = wave;
fprintf('Wind %g m/s, direction %g deg: solving RIS six-DOF motion ...\n', ...
    cfg.ris.windSpeed,cfg.ris.windAngleDeg);
m = RIS_float(cfg.ris);
step = round(cfg.comparison.frameStep/cfg.ris.dt);
assert(step>=1 && abs(step*cfg.ris.dt-cfg.comparison.frameStep)<1e-9, ...
    'wrd:FrameStep','comparison.frameStep must be a positive multiple of ris.dt.');
keys = unique([1:step:numel(m.t),numel(m.t)]);
rows = cfg.comparison.rotationXYDeg;
rx.nominal = cfg.receiver.center(:);
rx.position = repmat(rx.nominal,1,numel(m.t));
motions = cell(size(rows,1),1);
frames = cell(size(rows,1),numel(keys));
for row = 1:size(rows,1)
    cfg.geometry.rotationXYDeg = rows(row,:);
    motions{row} = wris_tx_motion(cfg,m.t);   % same TX model and seed as W_RIS
end
flatSea = wrd_flat_surface();
for frame = 1:numel(keys)
    k = keys(frame);
    sea = wrd_surface(m,m.t(k),cfg.surface,cfg.ris.center);
    [eta,slope] = RIS_wave_field(m,cfg.ris.center(:).',m.t(k));
    for row = 1:size(rows,1)
        cfg.geometry.rotationXYDeg = rows(row,:);
        g = rp_geometry(cfg,m,motions{row},rx,k);
        seed = cfg.mc.seed;   % common random numbers: same seed for all links/frames
        direct = wrd_direct_mc(cfg,g,sea,seed);
        flat = [];
        if cfg.comparison.flatBaseline
            flat = wrd_direct_mc(cfg,g,flatSea,seed);
        end
        ris = wrd_ris_mc(cfg,g,seed,sea);
        residual = [direct.energyResidual,ris.energyResidual];
        if ~isempty(flat), residual(end+1) = flat.energyResidual; end %#ok<AGROW>
        assert(max(abs(residual))<1e-9,'wrd:Energy', ...
            'MC energy ledger failed to close (U=%g, row %d, t=%g s).',cfg.ris.windSpeed,row,m.t(k));
        frames{row,frame} = struct('t',m.t(k),'absoluteTime',m.absoluteTime(k), ...
            'geometry',g,'etaM',eta,'slope',slope, ...
            'direct',direct,'flat',flat,'ris',ris, ...
            'bottomDryFraction',mean(m.bottomSeaGapM(:,k)>1e-4), ...
            'topWetFraction',mean(m.topSeaClearanceM(:,k)<0), ...
            'risTiltDeg',acosd(min(1,max(-1,m.normal(3,k)))));
    end
    if frame==1 || frame==numel(keys) || mod(frame,c.run.progressEvery)==0
        fprintf('  wind %g m/s | t = %5.2f s | frame %d/%d | elapsed %.1f min\n', ...
            cfg.ris.windSpeed,m.t(k),frame,numel(keys),toc(runTimer)/60);
    end
end
block = cell(size(rows,1),1);
for row = 1:size(rows,1)
    block{row} = struct('windMps',cfg.ris.windSpeed,'windDirectionDeg',cfg.ris.windAngleDeg, ...
        'txRotationXYDeg',rows(row,:), ...
        'nominalIncidenceDeg',acosd(cosd(rows(row,1))*cosd(rows(row,2))), ...
        'wave',m.wave,'significantWaveHeightM',4*sqrt(sum(abs(m.wave.coefficient).^2)/2), ...
        'risValidity',m.validity,'frames',vertcat(frames{row,:}));
end
end

function checkConfig(c)
validateattributes(c.mc.photons,{'numeric'},{'scalar','integer','>=',2});
validateattributes(c.mc.chunkSize,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.comparison.frameStep,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.comparison.conditions,{'numeric'},{'2d','ncols',2,'nonempty','real','finite'});
validateattributes(c.comparison.rotationXYDeg,{'numeric'},{'2d','ncols',2,'nonempty','real','finite'});
validateattributes(c.receiver.center,{'numeric'},{'vector','numel',3,'finite'});
validateattributes(c.receiver.normal,{'numeric'},{'vector','numel',3,'real','finite'});
validateattributes(c.receiver.radius,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.receiver.fovDeg,{'numeric'},{'scalar','positive','<=',90});
validateattributes(c.receiver.mapEdges,{'numeric'},{'vector','increasing'});
assert(norm(c.receiver.normal)>0,'wrd:Receiver','The RX normal must be nonzero.');
assert(c.receiver.center(3)>c.ris.thickness,'wrd:Receiver','Use a receiver above the RIS and the sea.');
validateattributes(c.surface.marchStep,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.surface.rootTolerance,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.surface.gridStep,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.surface.gridHalfWidth,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes(c.metrics.outageThresholdW,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.metrics.displayFloorW,{'numeric'},{'scalar','positive','finite'});
assert(islogical(c.comparison.flatBaseline) || isnumeric(c.comparison.flatBaseline), ...
    'wrd:Config','comparison.flatBaseline must be true or false.');
assert(c.ris.dt<=c.comparison.frameStep && c.comparison.frameStep<=c.ris.duration, ...
    'wrd:FrameStep','Need ris.dt <= comparison.frameStep <= ris.duration.');
for row = 1:size(c.comparison.rotationXYDeg,1)
    a = c.comparison.rotationXYDeg(row,:);
    assert(cosd(a(1))*cosd(a(2))>1e-6,'wrd:Direction', ...
        'TX axis [%g %g] deg must point upward toward the RIS lower face.',a);
end
end

function p = physicalConfig(c)
% Fields that change the transported numbers. Output/run/plot settings and the
% statistics thresholds (applied afterwards in wrd_summary) do not.
p = rmfield(c,intersect(fieldnames(c),{'output','run','plot','cache','root','profile','metrics'}));
end

function h = algorithmHash()
% Any edit to an algorithm file invalidates old checkpoints.
folder = fileparts(mfilename('fullpath'));
files = dir(fullfile(folder,'*.m'));
names = sort({files.name});
digest = java.security.MessageDigest.getInstance('SHA-256');
for j = 1:numel(names)
    digest.update(unicode2native(fileread(fullfile(folder,names{j})),'UTF-8'));
end
h = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
