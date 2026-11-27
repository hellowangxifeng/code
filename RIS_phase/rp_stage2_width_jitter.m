function result=rp_stage2_width_jitter(s)
%RP_STAGE2_WIDTH_JITTER Re-select actual feasible phase at EVERY jitter sigma.
c=s.config; d=s.distances{2}; samples=rp_link_samples(c);
candidateIds=find(d.summary.feasible); n=numel(candidateIds);
sigma=c.study.jitterGridM; rows=zeros(n*numel(sigma),9); row=0;
for j=candidateIds.'
    b=d.bank{j};
    for k=1:numel(sigma)
        opts=c.study; opts.positionSigmaM=sigma(k);
        h=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,samples,opts);
        m=rp_link_metrics(c,h,c.study.designPowerW); row=row+1;
        rows(row,:)=[j,k,sigma(k),100*b.spot.equivalentRadiusM,m.meanReceivedPowerW, ...
            m.aber,m.op,mean(vecnorm(h.positionOffsetM,2,2)>b.spot.r90M),d.summary.feasible(j)];
    end
end
T=array2table(rows,'VariableNames',{'candidate','sigmaId','sigmaM','actualRadiusCM', ...
    'meanReceivedPowerW','aber','op','probabilityOutsideR90','feasible'});
T.family=d.summary.family(T.candidate);
selection=zeros(numel(sigma)*5,9); count=0;
for k=1:numel(sigma)
    for objective=1:5
        pool=T(T.sigmaId==k & T.feasible>0,:);
        if objective<=2, pool=pool(pool.family=="robust",:); end
        if objective>=4, pool=pool(pool.family=="flat_top",:); end
        if objective==1 || objective==4
            [~,idx]=min(pool.aber);
        elseif objective==2 || objective==5
            [~,idx]=min(pool.op);
        else
            [~,idx]=max(pool.meanReceivedPowerW);
        end
        q=pool(idx,:); count=count+1;
        selection(count,:)=[sigma(k),objective,q.candidate,q.actualRadiusCM, ...
            q.meanReceivedPowerW,q.aber,q.op,q.probabilityOutsideR90,strcmp(q.family,"flat_top")];
    end
end
best=array2table(selection,'VariableNames',{'sigmaM','objectiveId','candidate','actualRadiusCM', ...
    'meanReceivedPowerW','aber','op','probabilityOutsideR90','selectedFlatTop'});
result=struct('table',T,'selection',best,'objectiveNames',{{'Robust ABER','Robust OP','Max mean power','Flat top: min ABER','Flat top: min OP'}}, ...
    'distanceM',d.distanceM,'designPowerW',c.study.designPowerW);
save(fullfile(c.output,'width_jitter.mat'),'result');
writetable(T,fullfile(c.output,'width_jitter.csv')); writetable(best,fullfile(c.output,'adaptive_widths.csv'));
end
