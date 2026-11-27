function key=rp_stage2_receiver_signature(s)
%RP_STAGE2_RECEIVER_SIGNATURE Bind scan caches to actual fields and receiver code.
c=s.config; cfg=rmfield(c,{'root','version','output','cache','plot'});
if isfield(cfg.study,'signature'), cfg.study=rmfield(cfg.study,'signature'); end
cfg.study=orderfields(cfg.study); cfg=orderfields(cfg);
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native([version,jsonencode(cfg)],'UTF-8'));
for j=1:numel(s.distances)
    d=s.distances{j}; digest.update(unicode2native(jsonencode(d.geometry),'UTF-8'));
    digest.update(typecast(double(d.comparisonIndices),'uint8'));
    for k=1:numel(d.bank)
        b=d.bank{k};
        digest.update(typecast(double(b.wave.irradianceWm2(:)),'uint8'));
        digest.update(unicode2native(jsonencode(struct('x',b.wave.x,'y',b.wave.y, ...
            'spot',b.spot,'family',b.family,'targetWidth',b.targetWidthM)),'UTF-8'));
    end
end
for name={'rp_stage2_sensitivity','rp_stage2_width_jitter','rp_stage2_turbulence_scan', ...
        'rp_stage2_select','rp_aperture_channel','rp_link_samples','rp_link_metrics', ...
        'rp_fading_response','rp_turbulence','rp_direct_wave','rp_direct_density','rp_direct_channel'}
    digest.update(unicode2native(fileread(which(name{1})),'UTF-8'));
end
key=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
