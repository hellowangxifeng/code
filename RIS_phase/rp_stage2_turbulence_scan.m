function result=rp_stage2_turbulence_scan(s)
%RP_STAGE2_TURBULENCE_SCAN Re-optimize at 1 W for each EGG/GG condition.
c=s.config; d=s.distances{2}; direct=rp_direct_wave(c,d.geometry);
names={'No RIS direct','Steering-only','Flat top','Robust ABER','Robust OP','Max mean power'};
rows=[]; choices=[]; kinds={'weak','strong'};
for ocean=1:2
    for atmosphere=1:2
        cc=c; cc.study.turbulence=kinds{ocean}; cc.study.atmosphere=kinds{atmosphere}; dc=d;
        for j=1:numel(dc.bank)
            dc.bank{j}.designMetrics=rp_link_metrics(cc,dc.bank{j}.channel,c.study.designPowerW);
        end
        dc=rp_stage2_select(cc,dc); indices=[2,dc.comparisonIndices(3:6)];
        bank=[{direct},dc.bank(indices)]; candidateIds=[0 indices];
        for objective=1:3
            id=dc.objectiveIndices(objective); b=dc.bank{id}; m=b.designMetrics;
            choices=[choices;ocean,atmosphere,objective,id,b.targetWidthM, ...
                100*b.spot.equivalentRadiusM,m.meanReceivedPowerW,m.aber,m.op]; %#ok<AGROW>
        end
        for scheme=1:6
            b=bank{scheme}; m=rp_link_metrics(cc,b.channel,10.^(c.study.powerDBm/10)/1000); n=numel(m.powerW);
            rows=[rows;repmat([ocean atmosphere scheme],n,1),c.study.powerDBm(:), ...
                m.meanReceivedPowerW(:),m.aber(:),m.op(:), ...
                repmat([candidateIds(scheme),100*b.spot.equivalentRadiusM,b.targetWidthM],n,1)]; %#ok<AGROW>
        end
    end
end
selection=array2table(choices,'VariableNames',{'oceanId','atmosphereId','objectiveId','candidate', ...
    'targetWidthM','actualRadiusCM','meanReceivedPowerW','aber','op'});
result=struct('table',array2table(rows,'VariableNames',{'oceanId','atmosphereId','schemeId', ...
    'transmitDBm','meanReceivedPowerW','aber','op','candidate','actualRadiusCM','targetWidthM'}), ...
    'selection',selection,'schemeNames',{names},'distanceM',7.5,'direct',direct, ...
    'note','Each EGG/GG condition reselects the phase at 1 W; that phase is frozen during its power sweep.');
save(fullfile(c.output,'turbulence_scan.mat'),'result','-v7.3');
writetable(result.table,fullfile(c.output,'turbulence_scan.csv'));
writetable(selection,fullfile(c.output,'turbulence_optimized_phases.csv'));
end
