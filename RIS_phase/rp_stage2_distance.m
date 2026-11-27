function d=rp_stage2_distance(c,z)
%RP_STAGE2_DISTANCE Propagate nominal families before communication selection.
rp_setup(); [g,mc,a]=rp_stage2_scene(c,z);
folder=fullfile(c.output,sprintf('z_%04.1f',z)); if ~isfolder(folder), mkdir(folder); end
diary(fullfile(folder,'run.log')); clean=onCleanup(@() diary('off'));
bank={rp_stage2_wave(c,g,mc,a,'unmodulated',0,0), ...
    rp_stage2_wave(c,g,mc,a,'steering',0,0)};
samples=rp_link_samples(c); baseline=rp_phase(c,g,a); screens=[];
for k=1:numel(c.study.families)
    family=c.study.families{k};
    for width=c.study.widthsM
        score=inf(size(c.study.centerMargins)); angles=score;
        for j=1:numel(c.study.centerMargins)
            margin=c.study.centerMargins(j);
            p=rp_finite_phase(c,g,a,family,width,margin,baseline);
            proxy=rp_proxy_wave(c,g,a,p); f=rp_spot_stats(proxy);
            h=rp_aperture_channel(c,g,proxy,f.centerXY,samples);
            m=rp_link_metrics(c,h,c.study.designPowerW,24);
            angles(j)=atan2d(norm(f.centerXY-g.pRis(1:2).'),z); score(j)=m.aber;
            screens=[screens;k,width,margin,angles(j),f.equivalentRadiusM,m.aber,m.op]; %#ok<AGROW>
        end
        score(angles>min(angles)+c.study.directionSlackDeg)=Inf;
        [~,best]=min(score);
        bank{end+1}=rp_stage2_wave(c,g,mc,a,family,width,c.study.centerMargins(best)); %#ok<AGROW>
        progress=struct('stage','nominal_waves','distanceM',z,'completed',numel(bank),'total',29);
        fid=fopen(fullfile(folder,'progress.json'),'w'); fprintf(fid,'%s',jsonencode(progress)); fclose(fid);
    end
end
d=struct('distanceM',z,'geometry',g,'mc',mc,'illumination',a,'bank',{bank},'screening',screens);
d=rp_stage2_select(c,d); save(fullfile(folder,'distance.mat'),'d','-v7.3');
end
