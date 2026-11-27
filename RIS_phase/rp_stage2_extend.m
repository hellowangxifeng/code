function d=rp_stage2_extend(c,d)
%RP_STAGE2_EXTEND Resolve width-boundary optima, with wider bank for sigma=120 mm.
rp_setup(); [g,mc,a]=rp_stage2_scene(c,d.distanceM);
widths=c.study.extendedWidthsM;
if abs(d.distanceM-7.5)<1e-9, widths=c.study.jitterExtendedWidthsM; end
if isfield(c.study,'lowerWidthsM'), widths=[c.study.lowerWidthsM,widths]; end
samples=rp_link_samples(c); baseline=rp_phase(c,g,a);
for family={'flat_top','robust'}
    for width=widths
        found=cellfun(@(b) strcmp(b.family,family{1}) && abs(b.targetWidthM-width)<1e-9,d.bank);
        if any(found), continue; end
        scores=inf(size(c.study.centerMargins)); angles=scores;
        for j=1:numel(c.study.centerMargins)
            p=rp_finite_phase(c,g,a,family{1},width,c.study.centerMargins(j),baseline);
            proxy=rp_proxy_wave(c,g,a,p); spot=rp_spot_stats(proxy);
            h=rp_aperture_channel(c,g,proxy,spot.centerXY,samples);
            m=rp_link_metrics(c,h,c.study.designPowerW,24);
            angles(j)=atan2d(norm(spot.centerXY-g.pRis(1:2).'),d.distanceM); scores(j)=m.aber;
        end
        scores(angles>min(angles)+c.study.directionSlackDeg)=Inf;
        [~,chosen]=min(scores);
        d.bank{end+1}=rp_stage2_wave(c,g,mc,a,family{1},width,c.study.centerMargins(chosen));
        rp_study_progress(c,sprintf('extended_widths_z%g',d.distanceM),numel(d.bank),29+2*numel(widths));
    end
end
d.bank=rp_stage2_sort_bank(d.bank); d=rp_stage2_select(c,d);
save(fullfile(c.output,sprintf('z_%04.1f',d.distanceM),'distance_extended.mat'),'d','-v7.3');
local=c; local.output=fullfile(c.output,sprintf('z_%04.1f',d.distanceM));
rp_study_progress(local,'width_bank_complete',numel(d.bank),numel(d.bank),'complete');
end
