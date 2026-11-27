function [r,files]=rp_fig12_run(s)
%RP_FIG12_RUN Three-cm weak/weak adaptive study and paired numerical checks.
arguments
    s (1,1) struct
end
c=s.config; d=s.distances{2}; c.rx.apertureRadius=.03;
c.study.turbulence='weak'; c.study.atmosphere='weak';
c.study.fovHalfAngleRad=.240; c.study.angleSigmaRad=.025;
c.study.angleBiasRad=[.025 0]; c.study.designPowerW=1;
c.study.jitterGridM=linspace(0,.12,31);
c.study.samples=32768; c.study.validationSamples=32768;
c.plot.figureInches=[2.67 2.37]; c.plot.fontName='Times New Roman';
c.plot.tickFontSize=10; c.plot.resolution=600;
folder=fullfile(c.output,'fig12_revision'); if ~isfolder(folder), mkdir(folder); end
signature=inputSignature(c,d);
file=fullfile(folder,'Fig12_data.mat');
if isfile(file)
    cached=load(file,'r');
    if isfield(cached.r,'complete') && cached.r.complete && strcmp(cached.r.signature,signature)
        r=cached.r; files=rp_fig12_plot(r); return
    end
end
ids=find(d.summary.feasible).';
assert(all(d.summary.r50M(ids)>=.02),'rp:FiniteCore','A selected candidate violates the core constraint.');
samples=rp_link_samples(c,c.study.samples,14012026);
sigma=c.study.jitterGridM; values=zeros(numel(ids)*numel(sigma),8); count=0;
for id=ids
    b=d.bank{id};
    for k=1:numel(sigma)
        opts=c.study; opts.positionSigmaM=sigma(k);
        h=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,samples,opts);
        m=rp_fig12_metrics(c,h,1);
        count=count+1;
        values(count,:)=[id,k,sigma(k),100*b.spot.equivalentRadiusM, ...
            m.meanReceivedPowerW,m.aber,m.op,h.outsideMapApertureFraction];
    end
    writeProgress(folder,'adaptive_candidates',find(ids==id,1),numel(ids));
end
T=array2table(values,'VariableNames',{'candidate','sigmaId','sigmaM','actualRadiusCM', ...
    'meanReceivedPowerW','aber','op','outsideMapFraction'});
T.family=d.summary.family(T.candidate); best=rp_fig12_select(T);
r=struct('config',c,'folder',folder,'signature',signature,'candidateTable',T, ...
    'selection',best,'candidateMetadata',d.summary,'distanceM',d.distanceM);
save(file,'r');
r.stability=pairedChecks(c,d,best,folder);
save(file,'r');
r.environment=environmentChecks(c,d,ids,folder);
save(file,'r');
[r.frozen,r.power,r.frozenDesigns]=frozenChecks(c,d,best);
r.aperture=apertureChecks(c,d,best,folder);
r.complete=true; save(file,'r'); exportTables(r);
files=rp_fig12_plot(r); rp_fig12_report(r);
writeProgress(folder,'complete',1,1);
end

function key=inputSignature(c,d)
digest=java.security.MessageDigest.getInstance('SHA-256');
cfg=rmfield(c,{'output','cache','plot','version','root'});
if isfield(cfg.study,'signature'), cfg.study=rmfield(cfg.study,'signature'); end
digest.update(unicode2native(jsonencode(cfg),'UTF-8'));
digest.update(unicode2native(jsonencode(d.geometry),'UTF-8'));
digest.update(unicode2native(jsonencode(d.summary(:,{'candidate','family','feasible','r50M'})),'UTF-8'));
for id=find(d.summary.feasible).'
    b=d.bank{id}; digest.update(typecast(double(b.wave.irradianceWm2(:)),'uint8'));
    digest.update(unicode2native(jsonencode(struct('x',b.wave.x,'y',b.wave.y,'spot',b.spot)),'UTF-8'));
end
for name={'rp_fig12_run','rp_fig12_select','rp_fig12_metrics','rp_aperture_channel', ...
        'rp_link_samples','rp_link_metrics','rp_fading_response','rp_turbulence'}
    digest.update(unicode2native(fileread(which(name{1})),'UTF-8'));
end
key=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end

function result=pairedChecks(c,d,best,folder)
repeats=8; sigma=c.study.jitterGridM;
raw=zeros(numel(sigma)*5*repeats,6); count=0;
for rep=1:repeats
    samples=rp_link_samples(c,32768,14112026+rep*997);
    for k=1:numel(sigma)
        selected=best(best.sigmaId==k,:); opts=c.study; opts.positionSigmaM=sigma(k);
        for id=unique(selected.candidate).'
            b=d.bank{id}; h=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,samples,opts);
            m=rp_fig12_metrics(c,h,1);
            for objective=selected.objectiveId(selected.candidate==id).'
                count=count+1;raw(count,:)=[rep,k,objective,id,m.aber,m.op];
            end
        end
    end
    writeProgress(folder,'independent_paired_replicates',rep,repeats);
end
R=array2table(raw,'VariableNames',{'replicate','sigmaId','objectiveId','candidate','aber','op'});
rows=zeros(numel(sigma)*2,10); count=0;
for k=1:numel(sigma)
    for metric=1:2
        field='aber'; if metric==2, field='op'; end
        robust=sortrows(R(R.sigmaId==k & R.objectiveId==metric,:),'replicate');
        flat=sortrows(R(R.sigmaId==k & R.objectiveId==metric+3,:),'replicate');
        delta=flat.(field)-robust.(field); half=tinv(.975,repeats-1)*std(delta)/sqrt(repeats);
        count=count+1;
        rows(count,:)=[sigma(k),metric,mean(robust.(field)),mean(flat.(field)),mean(delta), ...
            half,mean(delta)-half,mean(delta)+half,sum(delta>0),sum(delta<0)];
    end
end
result=struct('raw',R,'summary',array2table(rows,'VariableNames',{'sigmaM','metricId', ...
    'robustMean','flatMean','flatMinusRobust','pairedHalf95','lower95','upper95', ...
    'robustWinsOf8','flatWinsOf8'}),'samplesPerReplicate',32768,'replicates',repeats);
end

function T=environmentChecks(c,d,ids,folder)
rows=zeros(12*5,11); count=0; condition=0;
for jitter=[.02 .08 .12]
    for angle=[.025 .060]
        for fov=[.240 .060]
            condition=condition+1; cc=c; cc.study.positionSigmaM=jitter;
            cc.study.angleSigmaRad=angle; cc.study.fovHalfAngleRad=fov;
            samples=rp_link_samples(cc,16384,15012026); metrics=zeros(numel(ids),8);
            for k=1:numel(ids)
                id=ids(k); b=d.bank{id};
                h=rp_aperture_channel(cc,d.geometry,b.wave,b.spot.centerXY,samples);
                m=rp_fig12_metrics(cc,h,1);
                metrics(k,:)=[id,1,jitter,100*b.spot.equivalentRadiusM, ...
                    m.meanReceivedPowerW,m.aber,m.op,h.fovOutageExact];
            end
            pool=array2table(metrics,'VariableNames',{'candidate','sigmaId','sigmaM','actualRadiusCM', ...
                'meanReceivedPowerW','aber','op','fovFloor'});
            pool.family=d.summary.family(pool.candidate); chosen=rp_fig12_select(pool);
            for objective=1:5
                q=chosen(chosen.objectiveId==objective,:);count=count+1;
                rows(count,:)=[condition,jitter,angle,fov,objective,q.candidate, ...
                    q.actualRadiusCM,q.aber,q.op,q.meanReceivedPowerW,q.fovFloor];
            end
            writeProgress(folder,'environment_reselection',condition,12);
        end
    end
end
T=array2table(rows,'VariableNames',{'conditionId','sigmaM','angleSigmaRad','fovHalfAngleRad', ...
    'objectiveId','candidate','actualRadiusCM','aber','op','meanReceivedPowerW','fovFloor'});
end

function [T,P,designs]=frozenChecks(c,d,best)
nominal=best(abs(best.sigmaM-.08)<1e-10,:);
flat=find(d.summary.family=="flat_top" & abs(d.summary.targetWidthM-.10)<1e-10,1);
designs=[nominal.candidate(nominal.objectiveId==1),nominal.candidate(nominal.objectiveId==3),flat];
samples=rp_link_samples(c,65536,16012026); sigma=c.study.jitterGridM;
rows=zeros(3*numel(sigma),7); count=0; powerDBm=-5:1:45;
powerRows=zeros(3*numel(powerDBm),5); powerCount=0;
for scheme=1:3
    b=d.bank{designs(scheme)};
    for k=1:numel(sigma)
        opts=c.study; opts.positionSigmaM=sigma(k);
        h=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,samples,opts);
        m=rp_fig12_metrics(c,h,1);count=count+1;
        rows(count,:)=[scheme,designs(scheme),sigma(k),m.aber,m.op,m.meanReceivedPowerW, ...
            mean(vecnorm(h.positionOffsetM,2,2)>b.spot.r90M)];
        if abs(sigma(k)-.12)<1e-10
            curve=rp_fig12_metrics(c,h,10.^(powerDBm/10)/1000);
            for n=1:numel(powerDBm)
                powerCount=powerCount+1;
                powerRows(powerCount,:)=[scheme,powerDBm(n),curve.aber(n),curve.op(n),curve.meanReceivedPowerW(n)];
            end
        end
    end
end
T=array2table(rows,'VariableNames',{'schemeId','candidate','sigmaM','aber','op', ...
    'meanReceivedPowerW','probabilityOutsideR90'});
P=array2table(powerRows,'VariableNames',{'schemeId','powerDBm','aber','op','meanReceivedPowerW'});
end

function T=apertureChecks(c,d,best,folder)
samples=rp_link_samples(c,32768,17012026); rows=zeros(height(best),8); count=0;
for k=1:numel(c.study.jitterGridM)
    opts=c.study;opts.positionSigmaM=c.study.jitterGridM(k);
    selected=best(best.sigmaId==k,:);
    for id=unique(selected.candidate).'
        b=d.bank{id}; h=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,samples,opts);
        m=rp_fig12_metrics(c,h,1);
        high=opts; high.apertureRadialOrder=8; high.apertureAzimuthCount=32;
        h=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,samples,high);
        refined=rp_fig12_metrics(c,h,1);
        for objective=selected.objectiveId(selected.candidate==id).'
            count=count+1;rows(count,:)=[opts.positionSigmaM,objective,id,m.aber,refined.aber, ...
                m.op,refined.op,abs(refined.meanReceivedPowerW/m.meanReceivedPowerW-1)];
        end
    end
    writeProgress(folder,'aperture_quadrature_check',k,numel(c.study.jitterGridM));
end
T=array2table(rows,'VariableNames',{'sigmaM','objectiveId','candidate','aberBase','aberRefined', ...
    'opBase','opRefined','meanPowerRelativeChange'});
end

function exportTables(r)
writetable(r.candidateTable,fullfile(r.folder,'Fig12_candidates.csv'));
writetable(r.selection,fullfile(r.folder,'Fig12_adaptive_widths.csv'));
writetable(r.stability.raw,fullfile(r.folder,'Fig12_paired_replicates.csv'));
writetable(r.stability.summary,fullfile(r.folder,'Fig12_stability.csv'));
writetable(r.environment,fullfile(r.folder,'Fig12_environment.csv'));
writetable(r.aperture,fullfile(r.folder,'Fig12_aperture_check.csv'));
writetable(r.frozen,fullfile(r.folder,'Fig13_frozen_jitter.csv'));
writetable(r.power,fullfile(r.folder,'Fig14_frozen_power.csv'));
end

function writeProgress(folder,stage,completed,total)
fid=fopen(fullfile(folder,'progress.json'),'w'); clean=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(struct('stage',stage,'completed',completed,'total',total)));
end
