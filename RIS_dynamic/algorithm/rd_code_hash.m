function h = rd_code_hash(group)
%RD_CODE_HASH SHA-256 of the source files that produce one kind of number.
%   rd_code_hash('field')   RIS designs, complex field and their MC check
%   rd_code_hash('direct')  no-RIS rough-surface Monte Carlo
%   rd_code_hash(names)     only the listed functions (resolved with WHICH)
% Post-processing files (metrics, selection, ranking, plots) are never
% hashed, so editing them never invalidates an expensive checkpoint.
common = {'rd_scene','rd_motion','rd_receiver','rd_capture','rd_compress','rd_map', ...
    'rd_transport','rp_geometry','rp_link_samples','rp_fresnel','wris_mc','wris_reference', ...
    'wris_intersect','RIS_float','RIS_wave_components','RIS_wave_field','RIS_flotation', ...
    'RIS_rotation','RIS_radiation','RIS_sidewall_load','get_default_RIS_params', ...
    'get_default_sea_params','Hen_Green','Table_abc'};
if iscell(group)
    names = group;
else
    switch group
        case 'field'
            names = [common,{'rd_candidate','rd_phase','rd_rx_direction','rd_wave', ...
                'rd_pixel_intensity','rd_trace','rd_ris_packets','rp_mc','rp_illuminate', ...
                'rp_bsf_ballistic','rp_incident_field','rp_phase','rp_vertical_command', ...
                'rp_carrier_landing','rp_paraxial','rp_command','rp_trace','rp_phase_profile', ...
                'rp_layer_path','rp_receiver_window','rp_fourier_window','rp_deposit_power', ...
                'rp_spot_stats','rp_aperture_channel','rp_gauss_legendre'}];
        case 'direct'
            names = [common,{'rd_condition','rd_direct_mc','rd_surface','rd_crossing', ...
                'rd_to_plane','rd_rx_direction','rp_carrier_landing'}];
        otherwise
            error('rd:HashGroup','Unknown hash group "%s".',group);
    end
end
paths = cellfun(@which,names(:),'UniformOutput',false);
digest = java.security.MessageDigest.getInstance('SHA-256');
for j = 1:numel(paths)
    if isempty(paths{j}) || ~isfile(paths{j})
        digest.update(unicode2native(sprintf('missing:%s',names{j}),'UTF-8'));
    else
        digest.update(unicode2native(fileread(paths{j}),'UTF-8'));
    end
end
h = lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
