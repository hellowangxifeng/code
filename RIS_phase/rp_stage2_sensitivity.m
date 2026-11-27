function result=rp_stage2_sensitivity(s)
%RP_STAGE2_SENSITIVITY z=5 m; reuse geometry for both oceans and FOV cutoffs.
c=s.config; d=s.distances{1}; samples=rp_link_samples(c);
parameters={'positionSigmaM','fovHalfAngleRad','angleSigmaRad'};
values={c.study.positionSweepM,c.study.fovSweepRad,c.study.angleSweepRad}; rows=[];
uniqueIds=unique(d.comparisonIndices,'stable'); kinds={'weak','strong'};
for bankId=uniqueIds
    b=d.bank{bankId}; schemes=find(d.comparisonIndices==bankId);
    for parameter=1:3
        for value=values{parameter}
            opts=c.study; opts.(parameters{parameter})=value;
            if parameter~=2, opts.fovHalfAngleRad=c.study.sensitivityFovHalfAngleRad; end
            if parameter==2
                h=b.channel; accepted=h.angleErrorRad<=value & cos(h.angleErrorRad)>0;
                h.gain=h.withoutFOV.*accepted; h.fovOutageSample=mean(~accepted);
                if opts.angleSigmaRad>0
                    h.fovOutageExact=ncx2cdf((value/opts.angleSigmaRad)^2,2, ...
                        sum(opts.angleBiasRad.^2)/opts.angleSigmaRad^2,'upper');
                else
                    h.fovOutageExact=double(norm(opts.angleBiasRad)>value);
                end
            else
                h=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,samples,opts);
            end
            for ocean=1:2
                cc=c; cc.study.turbulence=kinds{ocean}; cc.study.atmosphere='weak';
                m=rp_link_metrics(cc,h,c.study.designPowerW);
                for scheme=schemes
                    rows=[rows;ocean,scheme,parameter,value,m.meanReceivedPowerW,m.aber,m.op]; %#ok<AGROW>
                end
            end
        end
    end
end
T=array2table(rows,'VariableNames',{'oceanId','schemeId','parameterId','parameterValue','meanReceivedPowerW','aber','op'});
T=sortrows(T,{'oceanId','schemeId','parameterId','parameterValue'});
result=struct('table',T,'schemeNames',{d.comparisonNames},'oceanKinds',{{'weak','strong'}}, ...
    'atmosphere','weak','baselineFovHalfAngleRad',c.study.sensitivityFovHalfAngleRad);
save(fullfile(c.output,'sensitivity.mat'),'result'); writetable(T,fullfile(c.output,'sensitivity.csv'));
end
