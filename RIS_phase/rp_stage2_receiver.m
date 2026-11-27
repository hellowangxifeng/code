function extra=rp_stage2_receiver(s)
%RP_STAGE2_RECEIVER Complete signed sensitivity, adaptive width and fading scans.
c=s.config; file=fullfile(c.output,'receiver_extra.mat'); signature=rp_stage2_receiver_signature(s);
if isfile(file)
    cached=load(file,'extra','signature','complete');
    if isfield(cached,'complete') && cached.complete && strcmp(cached.signature,signature)
        extra=cached.extra; return
    end
end
extra=struct(); complete=false;
rp_study_progress(c,'receiver_sensitivity',0,3);
extra.sensitivity=rp_stage2_sensitivity(s); save(file,'extra','signature','complete','-v7.3');
rp_study_progress(c,'width_jitter_reselection',1,3);
extra.widthJitter=rp_stage2_width_jitter(s); save(file,'extra','signature','complete','-v7.3');
rp_study_progress(c,'ocean_atmosphere_scan',2,3);
extra.turbulence=rp_stage2_turbulence_scan(s); complete=true;
save(file,'extra','signature','complete','-v7.3');
rp_study_progress(c,'receiver_scans_complete',3,3);
end
