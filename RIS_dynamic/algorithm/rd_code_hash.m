function h = rd_code_hash(names)
%RD_CODE_HASH SHA-256 of the source files that produce the transported numbers.
%   rd_code_hash()        RIS_dynamic physics files + reused project functions
%   rd_code_hash(names)   only the listed functions (resolved with WHICH)
% Post-processing / run-control files are excluded, so editing them never
% invalidates the expensive Monte Carlo checkpoints.
if nargin==0
    here = fileparts(mfilename('fullpath'));
    files = dir(fullfile(here,'rd_*.m'));
    skip = {'rd_run.m','rd_evaluate.m','rd_check.m','rd_config.m','rd_setup.m', ...
        'rd_cached.m','rd_code_hash.m','rd_physical.m','rd_schemes.m','rd_entry_table.m'};
    local = sort(setdiff({files.name},skip));
    paths = fullfile(here,local);
    names = {'rp_stage2_scene','rp_study_scene','rp_mc','rp_illuminate','rp_bsf_ballistic', ...
        'rp_incident_field','rp_phase','rp_vertical_command','rp_carrier_landing','rp_paraxial', ...
        'rp_finite_phase','rp_zero_phase','rp_command','rp_trace','rp_phase_profile', ...
        'rp_fresnel','rp_layer_path','rp_geometry','rp_link_metrics','rp_fading_response', ...
        'rp_turbulence','rp_link_samples','rp_gauss_legendre','wris_mc','wris_reference', ...
        'wris_intersect','wris_metrics','RIS_float','RIS_wave_components','RIS_wave_field', ...
        'RIS_flotation','RIS_rotation','RIS_radiation','RIS_sidewall_load', ...
        'get_default_RIS_params','get_default_sea_params','Hen_Green','Table_abc'};
    paths = [paths(:);cellfun(@which,names(:),'UniformOutput',false)];
else
    paths = cellfun(@which,names(:),'UniformOutput',false);
end
digest = java.security.MessageDigest.getInstance('SHA-256');
for j = 1:numel(paths)
    if isempty(paths{j}) || ~isfile(paths{j})
        digest.update(unicode2native(sprintf('missing:%d',j),'UTF-8'));
    else
        digest.update(unicode2native(fileread(paths{j}),'UTF-8'));
    end
end
h = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
