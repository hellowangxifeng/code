function rd_check(c)
%RD_CHECK Fail early, with a readable message, before a long run.
try
    tooOld = isMATLABReleaseOlderThan("R2022a");
catch
    tooOld = true;
end
if tooOld
    error('rd:Release','RIS_dynamic needs MATLAB R2022a or newer (tiledlayout, clim, exportgraphics).');
end
if isempty(which('sobolset')) || isempty(which('ncx2cdf'))
    error('rd:Toolbox',['The Statistics and Machine Learning Toolbox is required ' ...
        '(sobolset, ncx2cdf), exactly as in RIS_phase.']);
end
needed = {'rp_stage2_config','rp_stage2_scene','rp_phase','rp_finite_phase','rp_trace', ...
    'rp_phase_profile','rp_link_metrics','rp_fading_response','rp_link_samples','rp_fresnel', ...
    'rp_geometry','wris_mc','wris_reference','wris_intersect','RIS_float', ...
    'RIS_wave_components','RIS_wave_field','RIS_flotation','get_default_sea_params', ...
    'Hen_Green','Table_abc'};
missing = needed(cellfun(@(n) isempty(which(n)),needed));
if ~isempty(missing)
    error('rd:Dependencies','Functions not found on the path: %s',strjoin(missing,', '));
end
assert(strcmp(c.phase.model,'linear'),'rd:Phase', ...
    'RIS_dynamic freezes affine LC phases; RIS_phase c.phase.model must be ''linear''.');
assert(c.source.power==1,'rd:Power','Use a 1 W reference source; sweep c.study.powerDBm instead.');
validateattributes(c.scene.distancesM,{'numeric'},{'vector','positive','finite','increasing'});
validateattributes(c.scene.incidenceDeg,{'numeric'},{'scalar','>',0,'<',80});
validateattributes(c.sea.windsMps,{'numeric'},{'vector','positive','finite','increasing'});
counts = [c.dynamic.photons c.dynamic.chunkSize c.dynamic.samples c.dynamic.calmFrames ...
    c.training.photons c.training.samples c.incidence.photons c.incidence.samples ...
    c.receiver.pilotFrames c.receiver.pilotPhotons c.study.samples c.maps.bins];
validateattributes(counts,{'numeric'},{'integer','positive'});
steps = [c.dynamic.frameStepS c.training.frameStepS c.incidence.frameStepS];
for s = steps
    n = round(s/c.ris.dt);
    assert(n>=1 && abs(n*c.ris.dt-s)<1e-9 && s<=c.ris.duration,'rd:FrameStep', ...
        'Frame steps must be multiples of c.ris.dt and not longer than c.ris.duration.');
end
validateattributes([c.surface.gridStep c.surface.marchStep c.surface.rootTolerance], ...
    {'numeric'},{'positive','finite'});
validateattributes(c.surface.gridHalfWidth,{'numeric'},{'scalar','nonnegative','finite'});
validateattributes([c.compress.positionM c.compress.directionCos c.capture.maxPairs], ...
    {'numeric'},{'positive','finite'});
assert(ismember(lower(c.training.objective),{'aber','op'}),'rd:Training', ...
    'c.training.objective must be ''aber'' or ''op''.');
assert(size(c.metrics.turbulence,2)==2 && iscellstr(c.metrics.turbulence),'rd:Turbulence', ...
    'c.metrics.turbulence must be an n-by-2 cell of {ocean, atmosphere} names.');
kinds = {'none','weak','strong'};
assert(all(ismember(c.metrics.turbulence(:),kinds)),'rd:Turbulence','Use none, weak or strong.');
if c.incidence.enabled
    validateattributes(c.incidence.anglesDeg,{'numeric'},{'vector','>',0,'<',80,'increasing'});
    assert(c.incidence.windMps>0,'rd:Incidence','c.incidence.windMps must be positive.');
    if c.training.enabled && ~ismember(c.incidence.windMps,c.sea.windsMps)
        warning('rd:Incidence',['c.incidence.windMps is not in c.sea.windsMps; the incidence ' ...
            'sweep will use Robust ABER in place of the sea-aware design.']);
    end
end
end
