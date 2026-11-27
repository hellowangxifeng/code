function results = run_wris_dynamic(input)
%RUN_WRIS_DYNAMIC Calculate and save paired rough-sea and moving-RIS data.
%   run_wris_dynamic() uses the standard settings.
%   run_wris_dynamic('smoke') runs a small end-to-end case.
%   run_wris_dynamic(c) accepts an edited WRD_CONFIG structure.
if nargin==0, input = 'production'; end
if isstruct(input)
    wrd_config('smoke'); % register the read-only project dependencies
    c = input;
else
    c = wrd_config(input);
end
validateattributes(c.mc.photons,{'numeric'},{'scalar','integer','>=',2});
validateattributes(c.comparison.frameStep,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.receiver.center,{'numeric'},{'vector','numel',3,'finite'});
validateattributes(c.receiver.normal,{'numeric'},{'vector','numel',3,'real','finite'});
assert(norm(c.receiver.normal)>0,'wrd:Receiver','The RX normal must be nonzero.');
assert(c.receiver.center(3)>0,'wrd:Receiver','Use an above-sea receiver.');
validateattributes(c.surface.marchStep,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.surface.rootTolerance,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.surface.gridStep,{'numeric'},{'scalar','positive','finite'});
validateattributes(c.surface.gridHalfWidth,{'numeric'},{'scalar','nonnegative','finite'});

rows = c.comparison.rotationXYDeg;
conditions = c.comparison.conditions;
cases = cell(size(rows,1),size(conditions,1));
seaDefault = get_default_sea_params();
for seaId = 1:size(conditions,1)
    cfg = c;
    cfg.ris.windSpeed = conditions(seaId,1);
    cfg.ris.windAngleDeg = conditions(seaId,2);
    wave = RIS_wave_components(cfg.ris,cfg.sea.quadrature);
    wave.coefficient = wave.coefficient*sqrt(cfg.sea.neumannConstant/seaDefault.neumannConstant);
    cfg.ris.wave = wave;
    m = RIS_float(cfg.ris);
    step = round(cfg.comparison.frameStep/cfg.ris.dt);
    assert(step>=1 && abs(step*cfg.ris.dt-cfg.comparison.frameStep)<1e-9, ...
        'wrd:FrameStep','Use a frame step that is a positive multiple of the motion output step.');
    [motions,frames] = deal(cell(size(rows,1),1));
    keys = unique([1:step:numel(m.t),numel(m.t)]);
    rx.nominal = cfg.receiver.center(:);
    rx.position = repmat(rx.nominal,1,numel(m.t));
    for row = 1:size(rows,1)
        cfg.geometry.rotationXYDeg = rows(row,:);
        motions{row} = wris_tx_motion(cfg,m.t);
        frames{row} = cell(numel(keys),1);
    end
    for frame = 1:numel(keys)
        k = keys(frame);
        sea = wrd_surface(m,m.t(k),cfg.surface,cfg.ris.center);
        for row = 1:size(rows,1)
            cfg.geometry.rotationXYDeg = rows(row,:);
            g = rp_geometry(cfg,m,motions{row},rx,k);
            % Share the random seed between links and between frames.
            seed = cfg.mc.seed;
            direct = wrd_direct_mc(cfg,g,sea,seed);
            ris = wrd_ris_mc(cfg,g,seed);
            assert(max(abs([direct.energyResidual,ris.energyResidual]))<1e-9, ...
                'wrd:Energy','The MC transport energy ledger failed to close.');
            [eta,slope] = RIS_wave_field(m,cfg.ris.center(:).',m.t(k));
            frames{row}{frame} = struct('t',m.t(k),'absoluteTime',m.absoluteTime(k), ...
                'geometry',g,'etaM',eta,'slope',slope,'direct',direct,'ris',ris, ...
                'bottomDryFraction',mean(m.bottomSeaGapM(:,k)>1e-4), ...
                'topWetFraction',mean(m.topSeaClearanceM(:,k)<0));
        end
        if frame==1 || frame==numel(keys) || mod(frame,5)==0
            fprintf('Wind %.1f m/s, time %.2f s, frame %d/%d\n', ...
                cfg.ris.windSpeed,m.t(k),frame,numel(keys));
        end
    end
    for row = 1:size(rows,1)
        cases{row,seaId} = struct('windMps',cfg.ris.windSpeed, ...
            'windDirectionDeg',cfg.ris.windAngleDeg,'wave',m.wave, ...
            'txRotationXYDeg',rows(row,:),'frames',vertcat(frames{row}{:}));
    end
end
results = struct('config',c,'cases',{cases},'summary',buildSummary(cases));
if ~isfolder(c.output), mkdir(c.output); end
save(fullfile(c.output,'wris_dynamic.mat'),'results','-v7');
writetable(results.summary,fullfile(c.output,'wris_dynamic_summary.csv'));
fprintf('Saved paired %d cases to %s\n',numel(cases),c.output);
end

function tbl = buildSummary(cases)
rows = cell(numel(cases),1);
for id = 1:numel(cases)
    item = cases{id};
    f = item.frames;
    direct = [f.direct];
    ris = [f.ris];
    time = [f.t].';
    rows{id} = table(repmat(item.windMps,numel(f),1), ...
        repmat(item.windDirectionDeg,numel(f),1), ...
        repmat(item.txRotationXYDeg(1),numel(f),1),repmat(item.txRotationXYDeg(2),numel(f),1), ...
        time,[f.absoluteTime].',reshape([f.etaM],[],1),[f.bottomDryFraction].',[f.topWetFraction].', ...
        [direct.rxPowerW].',[direct.powerSEW].',[direct.hitCount].', ...
        [direct.tirPowerW].',[direct.incidentPowerW].',[direct.transmittedPowerW].', ...
        [direct.ballisticPowerW].',[direct.scatterPowerW].',[direct.waveBlockedPowerW].', ...
        [ris.rxPowerW].',[ris.powerSEW].',[ris.hitCount].',[ris.incidentPowerW].', ...
        [ris.activePowerW].',[ris.transmittedPowerW].',[ris.ballisticPowerW].',[ris.scatterPowerW].', ...
        [ris.saturatedPowerFraction].',[direct.energyResidual].',[ris.energyResidual].', ...
        'VariableNames',{'windMps','windDirectionDeg','txRotationXDeg','txRotationYDeg', ...
        'tS','absoluteTimeS','etaM','risBottomDrySampleFraction','risTopWetSampleFraction', ...
        'directRxW','directRxSEW','directRxHits','directTirW','directSeaIncidentW','directTransmittedW', ...
        'directRxBallisticW','directRxScatterW','directWaveBlockedW','risRxW','risRxSEW', ...
        'risRxHits','risBoardIncidentW','risActiveW','risTransmittedW','risRxBallisticW','risRxScatterW', ...
        'risSaturatedFraction','directEnergyResidual','risWaterEnergyResidual'});
end
tbl = vertcat(rows{:});
end
