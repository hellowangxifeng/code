function v=rp_stage2_validate(s,extra)
%RP_STAGE2_VALIDATE Physics, EGG/GG convergence, independent jitter and selections.
c=s.config; samples=rp_link_samples(c,c.study.validationSamples,c.study.seed+703);
receiverCache=load(fullfile(c.output,'receiver_extra.mat'),'signature','complete');
assert(receiverCache.complete && strcmp(receiverCache.signature,rp_stage2_receiver_signature(s)), ...
    'rp:ReceiverProvenance','Receiver scans are not bound to current fields and parameters.');
energy=[]; statistics=[];
for j=1:3
    d=s.distances{j};
    assert(all(d.summary.feasible(d.objectiveIndices)),'rp:Feasibility','Selected phase is infeasible.');
    assert(all(cellfun(@(b) b.spot.r50M>=.02,d.bank(d.objectiveIndices))), ...
        'rp:PointFocus','Half the energy is concentrated within a 2 cm disk.');
    for k=1:numel(d.bank)
        b=d.bank{k}; w=b.wave; p=b.phase;
        limit=floor(c.device.maxAngleDeg/c.device.angleStepDeg)*c.device.angleStepDeg;
        assert(max(abs(p.commandDeg),[],'all')<=limit+1e-12,'rp:Command','Command limit violated.');
        assert(max(abs(p.commandDeg/.021-round(p.commandDeg/.021)),[],'all')<1e-9,'rp:Quantization','Angle is not quantized.');
        assert(max(abs(p.gradientOverK0-c.optics.nAir*sind(p.commandDeg)),[],'all')<1e-12,'rp:Gradient','Gradient mismatch.');
        assert(all(p.phaseAtReferenceRad>=0 & p.phaseAtReferenceRad<2*pi),'rp:PhaseWrap','Phase outside 0--2pi.');
        assert(all(isfinite(w.irradianceWm2) & w.irradianceWm2>=0,'all'),'rp:Irradiance','Invalid intensity.');
        residual=w.pupilPowerW+w.interfaceLossW+w.omittedBallisticPowerW-d.illumination.activeDirectPowerW;
        excess=(sum(w.ballisticPixelMeanWm2,'all')*w.cellAreaM2-w.pupilPowerW)/max(w.pupilPowerW,realmin);
        assert(abs(residual)<1e-10 && excess<1e-3,'rp:Energy','Optical energy tolerance exceeded.');
        assert(w.diagonalMapPowerW<=w.pupilPowerW+1e-9,'rp:Parseval','Diagonal energy exceeds pupil.');
        assert(w.unresolvedQuadraturePowerFraction<c.detector.maximumUnresolvedPowerFraction,'rp:Quadrature','Unresolved optical power.');
        curve=rp_link_metrics(c,b.channel,10.^(c.study.powerDBm/10)/1000);
        assert(all(diff(curve.aber)<=1e-12) && all(diff(curve.op)<=1e-12),'rp:Monotonicity','BER or OP increases with power.');
        energy=[energy;j,k,residual,excess,w.diagonalMapPowerW/max(w.pupilPowerW,realmin),w.scatteredMapPowerW,w.unresolvedQuadraturePowerFraction]; %#ok<AGROW>
    end
    for k=unique(d.comparisonIndices)
        b=d.bank{k}; h=rp_aperture_channel(c,d.geometry,b.wave,b.spot.centerXY,samples);
        independent=rp_link_metrics(c,h,c.study.designPowerW,64);
        low=rp_link_metrics(c,b.channel,c.study.designPowerW,32);
        high=rp_link_metrics(c,b.channel,c.study.designPowerW,64);
        statistics=[statistics;j,k,abs(independent.meanReceivedPowerW/low.meanReceivedPowerW-1), ...
            abs(independent.aber-low.aber),abs(independent.op-low.op), ...
            abs(high.aber-low.aber),abs(high.op-low.op),h.outsideMapApertureFraction]; %#ok<AGROW>
    end
end
energy=array2table(energy,'VariableNames',{'distanceId','candidate','balanceW','coherentExcess', ...
    'diagonalCapture','scatterMapW','unresolvedFraction'});
statistics=array2table(statistics,'VariableNames',{'distanceId','candidate','relativeMeanPowerChange', ...
    'independentABERChange','independentOPChange','orderABERChange','orderOPChange','outsideMapApertureFraction'});
assert(max(statistics.relativeMeanPowerChange)<.02 && max(statistics.independentABERChange)<.002 && ...
    max(statistics.independentOPChange)<.005 && max(statistics.orderABERChange)<.001 && ...
    max(statistics.orderOPChange)<.002,'rp:Statistics','Independent EGG/GG statistics need refinement.');
assert(all(cellfun(@(q) q.passed,extra.optical)),'rp:Optical','Optical refinement failed.');
powerCheck=extra.optical{5};
assert(powerCheck.refinedSpot.r50M>=.02 && powerCheck.refinedSpot.minorRadiusM>=c.study.minimumWidthM, ...
    'rp:PowerCoreResolution','Power-optimal core violates minimum width after refinement.');
assert(extra.optical{6}.refinedSpot.r50M>=.02,'rp:AdaptiveCoreResolution', ...
    'Low-jitter robust core violates the 2 cm half-energy radius after refinement.');
% Independent closed-form moments distinguish atmospheric Gamma-Gamma from
% the generalized-gamma component inside ocean EGG.
moments=[]; kinds={'weak','strong'};
for ocean=1:2
    for atmosphere=1:2
        cc=c; cc.study.turbulence=kinds{ocean}; cc.study.atmosphere=kinds{atmosphere};
        t=rp_turbulence(cc,64); shapes=t.ggShapes;
        expected=1+1/shapes(1)+1/shapes(2)+1/prod(shapes);
        meanA=sum(t.atmosphereValues.*t.atmosphereWeights);
        secondA=sum(t.atmosphereValues.^2.*t.atmosphereWeights);
        assert(abs(meanA-1)<1e-11 && abs(secondA-expected)<1e-10,'rp:GGMoments','GG moments fail.');
        response=rp_fading_response(cc,64); meanJoint=response.distributionMean;
        assert(abs(meanJoint-1)<5e-6 && abs(response.distributionMass-1)<1e-10, ...
            'rp:EGGMoments','Tail-resolving EGG/GG distribution mass or mean inaccurate.');
        moments=[moments;ocean,atmosphere,shapes,meanA,secondA,expected,meanJoint]; %#ok<AGROW>
    end
end
moments=array2table(moments,'VariableNames',{'oceanId','atmosphereId','alpha','beta','GGMean','GGSecondMoment','exactGGSecondMoment','jointMeanLogGrid'});
fadingCheck=readtable(fullfile(c.root,'audit','stage2','joint_fading_verification.csv'));
assert(max(fadingCheck.ber64QMCAbsoluteError)<1e-5 && max(fadingCheck.op64QMCAbsoluteError)<5e-4, ...
    'rp:IndependentFading','Independent inverse-CDF fading verification failed.');
T=extra.widthJitter.table; best=extra.widthJitter.selection;
for j=1:height(best)
    q=best(j,:); pool=T(abs(T.sigmaM-q.sigmaM)<1e-12 & T.feasible>0,:);
    if q.objectiveId<=2, pool=pool(pool.family=="robust",:); end
    if q.objectiveId>=4, pool=pool(pool.family=="flat_top",:); end
    if q.objectiveId==1 || q.objectiveId==4
        assert(abs(q.aber-min(pool.aber))<1e-12,'rp:WidthChoice','Incorrect per-sigma ABER selection.');
    elseif q.objectiveId==2 || q.objectiveId==5
        assert(abs(q.op-min(pool.op))<1e-12,'rp:WidthChoice','Incorrect per-sigma OP selection.');
    else
        assert(abs(q.meanReceivedPowerW-max(pool.meanReceivedPowerW))<1e-12,'rp:WidthChoice','Incorrect per-sigma power selection.');
    end
end
for j=1:2
    q=extra.angles{j}.table;
    direct=q(q.schemeId==1 & q.incidenceDeg>=60,:);
    assert(all(direct.meanReceivedPowerW<1e-8 & direct.op>1-1e-8 & direct.aber>.5-1e-8), ...
        'rp:TIR','Direct-path TIR behavior is missing.');
    steered=q(q.schemeId==2 & q.incidenceDeg==60,:);
    assert(steered.meanReceivedPowerW>1e-6,'rp:RISRecovery','RIS has no propagated output above direct critical angle.');
    maps=extra.angles{j}.maps;
    for angle=1:size(maps,1)
        for scheme=2:6
            b=maps{angle,scheme}; p=b.phase; w=b.wave;
            assert(max(abs(p.commandDeg),[],'all')<=15.414+1e-12 && ...
                max(abs(p.commandDeg/.021-round(p.commandDeg/.021)),[],'all')<1e-9, ...
                'rp:AnglePhase','Invalid realized phase in incidence scan.');
            assert(all(p.phaseAtReferenceRad>=0 & p.phaseAtReferenceRad<2*pi),'rp:AngleWrap','Angle-scan piston is not wrapped.');
            balance=w.pupilPowerW+w.interfaceLossW+w.omittedBallisticPowerW-w.inputDirectPowerW;
            assert(abs(balance)<1e-10 && w.diagonalMapPowerW<=w.pupilPowerW+1e-9, ...
                'rp:AngleEnergy','Incidence-scan energy accounting failed.');
            assert(w.unresolvedQuadraturePowerFraction<c.detector.maximumUnresolvedPowerFraction, ...
                'rp:AngleQuadrature','Unresolved incidence-scan optical power.');
        end
    end
end
direct=extra.turbulence.direct;
directCapture=direct.spot.powerW/max(direct.wave.angularCubaturePowerW,realmin);
assert(abs(directCapture-1)<.02,'rp:DirectEnergy','Continuous direct-beam Jacobian/cubature disagree.');
adaptive=rp_stage2_validate_adaptive(s,extra);
v=struct('status','numerical_passed_visual_review_pending','modelSignature',c.study.signature, ...
    'statistics',statistics,'energy',energy,'fadingMoments',moments,'directMapCapture',directCapture, ...
    'adaptive',adaptive, ...
    'optical', {cellfun(@(q) rmfield(q,'refinedWave'),extra.optical,'UniformOutput',false)});
save(fullfile(c.output,'validation.mat'),'v');
writetable(statistics,fullfile(c.output,'validation_statistics.csv')); writetable(energy,fullfile(c.output,'validation_energy.csv'));
writetable(moments,fullfile(c.output,'validation_fading.csv'));
fid=fopen(fullfile(c.output,'validation.json'),'w'); clean=onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(v,PrettyPrint=true));
end
