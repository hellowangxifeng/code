function result=rp_stage2_validate_adaptive(s,extra)
%RP_STAGE2_VALIDATE_ADAPTIVE Independently evaluate EVERY reselected sigma point.
c=s.config; d=s.distances{2}; B=extra.widthJitter.selection;
samples=rp_link_samples(c,c.study.validationSamples,c.study.seed+703); rows=[];
for sigma=unique(B.sigmaM).'
    q=B(abs(B.sigmaM-sigma)<1e-12,:); opts=c.study; opts.positionSigmaM=sigma;
    for id=unique(q.candidate).'
        b=d.bank{id}; h=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,samples,opts);
        m=rp_link_metrics(c,h,c.study.designPowerW,64);
        for idx=find(q.candidate==id).'
            v=q(idx,:);
            rows=[rows;sigma,v.objectiveId,id,v.actualRadiusCM,v.aber,m.aber, ...
                v.op,m.op,abs(m.meanReceivedPowerW/v.meanReceivedPowerW-1)]; %#ok<AGROW>
        end
    end
end
T=array2table(rows,'VariableNames',{'sigmaM','objectiveId','candidate','actualRadiusCM', ...
    'aberDesign','aberIndependent','opDesign','opIndependent','meanPowerRelativeDifference'});
T=sortrows(T,{'sigmaM','objectiveId'});
assert(max(abs(T.aberDesign-T.aberIndependent))<.002 && ...
    max(abs(T.opDesign-T.opIndependent))<.005 && max(T.meanPowerRelativeDifference)<.02, ...
    'rp:AdaptiveConvergence','Adaptive curves need more position/angular samples.');
comparison=[];
for sigma=unique(T.sigmaM).'
    q=T(abs(T.sigmaM-sigma)<1e-12,:);
    r=q(q.objectiveId==1,:); o=q(q.objectiveId==2,:);
    f=q(q.objectiveId==4,:); g=q(q.objectiveId==5,:);
    comparison=[comparison;sigma,r.aberIndependent,f.aberIndependent, ...
        o.opIndependent,g.opIndependent,f.aberIndependent-r.aberIndependent, ...
        g.opIndependent-o.opIndependent,f.aberIndependent<r.aberIndependent-1e-6, ...
        g.opIndependent<o.opIndependent-1e-5]; %#ok<AGROW>
end
comparison=array2table(comparison,'VariableNames',{'sigmaM','robustABER','flatBestABER', ...
    'robustOP','flatBestOP','flatMinusRobustABER','flatMinusRobustOP','flatLowerABER','flatLowerOP'});
result=struct('table',T,'comparison',comparison,'samples',c.study.validationSamples, ...
    'maxABERChange',max(abs(T.aberDesign-T.aberIndependent)), ...
    'maxOPChange',max(abs(T.opDesign-T.opIndependent)));
writetable(T,fullfile(c.output,'validation_adaptive.csv'));
writetable(comparison,fullfile(c.output,'flat_top_advantage.csv'));
save(fullfile(c.output,'validation_adaptive.mat'),'result');
end
