function data=rp_fixed_receiver_collect(c,d,bank,unmodulated,policy)
%RP_FIXED_RECEIVER_COLLECT Fixed-RX performance and realized phase schedules.
assert(ismember(policy,{'per_power','fixed_design'}),'rp:FixedPolicy','Unknown phase policy.');
target=bank{1}.receiverCenterXY; powerDBm=0:.5:100;
powerW=10.^(powerDBm/10)/1000; n=numel(bank); nP=numel(powerW);
aber=zeros(n,nP); outage=aber; received=aber; rows=zeros(n,8);
direction=bank{1}.channel.nominalDirection;
for k=1:n
    b=bank{k}; h=b.channel;
    assert(norm(h.nominalCenterXY-target)<1e-12 && norm(h.nominalDirection-direction)<1e-12, ...
        'rp:FixedFrame','Schemes do not share a fixed receiver frame.');
    assert(h.outsideMapApertureFraction<1e-5,'rp:FixedMap','RX aperture leaves the computed field.');
    m=rp_fixed_receiver_metrics(c,h,powerW); aber(k,:)=m.aber; outage(k,:)=m.op; received(k,:)=m.meanReceivedPowerW;
    p=b.phase; limit=floor(c.device.maxAngleDeg/c.device.angleStepDeg)*c.device.angleStepDeg;
    assert(max(abs(p.commandDeg),[],'all')<=limit+1e-12 && ...
        max(abs(p.commandDeg/.021-round(p.commandDeg/.021)),[],'all')<1e-9,'rp:FixedPhase','Illegal fixed-target phase.');
    f=b.spot; rows(k,:)=[k,b.targetWidthM,f.equivalentRadiusM,f.minorRadiusM,f.axisRatio,f.r50M,f.r90M,b.wave.unresolvedQuadraturePowerFraction];
end
T=array2table(rows,'VariableNames',{'candidate','targetWidthM','actualRadiusM','minorRadiusM', ...
    'axisRatio','r50M','r90M','unresolvedFraction'});
T.family=string(cellfun(@(b)b.family,bank,'UniformOutput',false)).';
T.feasible=T.minorRadiusM>=c.study.minimumWidthM & T.actualRadiusM>=c.study.minimumWidthM & ...
    T.axisRatio<=c.study.maximumAxisRatio & T.r50M>=.02;
assert(all(T.unresolvedFraction<c.detector.maximumUnresolvedPowerFraction),'rp:FixedQuadrature','Unresolved propagation.');
feasible=T.feasible; feasible(1)=true; % steering is an explicit incumbent.
flat=feasible & T.family=="flat_top";
assert(any(flat),'rp:FixedFlat','No feasible flat-top candidate.');
designId=find(abs(powerW-c.study.designPowerW)==min(abs(powerW-c.study.designPowerW)),1);
choice=zeros(6,nP,2); % one phase per scheme/power; Flat is optimized per figure.
for metric=1:2
    choice(2,:,metric)=1;
    if metric==1, choice(3,:,metric)=choosePhase(aber,outage,received,flat);
    else, choice(3,:,metric)=choosePhase(outage,aber,received,flat); end
    choice(4,:,metric)=choosePhase(aber,outage,received,feasible);
    choice(5,:,metric)=choosePhase(outage,aber,received,feasible);
    score=received; score(~feasible,:)=-Inf; [~,choice(6,:,metric)]=max(score,[],1);
    if strcmp(policy,'fixed_design'), choice(:,:,metric)=repmat(choice(:,designId,metric),1,nP); end
end
curves=cell(1,6); optimistic=rp_fixed_receiver_metrics(c,unmodulated.channelUpper,powerW);
scatterOnly=rp_fixed_receiver_metrics(c,unmodulated.channelLower,powerW);
curves{1}=optimistic; exportRows=zeros(12*nP,11); count=0;
for scheme=1:6
    if scheme>1
        m=struct();
        for metric=1:2
            ids=choice(scheme,:,metric); linear=sub2ind([n nP],ids,1:nP);
            if metric==1, m.aber=aber(linear); m.meanReceivedPowerW=received(linear); else, m.op=outage(linear); end
        end
        curves{scheme}=m;
    end
    for metric=1:2
        for j=1:nP
            id=choice(scheme,j,metric); val=curves{scheme}.aber(j);
            if metric==2, val=curves{scheme}.op(j); end
            width=NaN; radius=NaN; if id>0, width=T.targetWidthM(id); radius=T.actualRadiusM(id); end
            count=count+1;
            exportRows(count,:)=[scheme,j,powerDBm(j),metric,id,width,radius,val,target,d.distanceM];
        end
    end
end
curveTable=array2table(exportRows,'VariableNames',{'schemeId','powerId','transmitDBm','metricId', ...
    'candidate','targetWidthM','actualRadiusM','metricValue','receiverXM','receiverYM','distanceM'});
folder=fullfile(c.output,'fixed_receiver');
writetable(T,fullfile(folder,'candidates.csv')); writetable(curveTable,fullfile(folder,'curves_and_phase_schedule.csv'));
feasibleABER=aber; feasibleABER(~feasible,:)=Inf;
feasibleOP=outage; feasibleOP(~feasible,:)=Inf;
if strcmp(policy,'per_power')
    assert(max(abs(curves{4}.aber-min(feasibleABER,[],1)))<1e-12);
    assert(max(abs(curves{5}.op-min(feasibleOP,[],1)))<1e-12);
end
data=struct('powerDBm',powerDBm,'curves',{curves},'schemeNames', ...
    {{'U-RIS optimistic envelope','Steering-only','Flat top','Robust ABER','Robust OP','Max mean power'}}, ...
    'displayDBm',[0 100],'figureSizeInches',[2.67 2.37],'fontPt',10,'fontName','Times New Roman', ...
    'fovHalfAngleRad',c.study.fovHalfAngleRad,'fovOutageExact',optimistic.fovOutageExact, ...
    'aberFloor',optimistic.fovOutageExact/2,'opFloor',optimistic.fovOutageExact, ...
    'distanceM',d.distanceM,'receiverCenterXY',target,'receiverDirection',direction, ...
    'positionSigmaM',c.study.positionSigmaM,'angleSigmaRad',c.study.angleSigmaRad, ...
    'phasePolicy',policy,'choice',choice,'candidateTable',T,'curveTable',curveTable, ...
    'sourceModelSignature',c.study.signature,'unmodulatedScatterOnly',scatterOnly, ...
    'unmodulatedTailEnvelopeWm2',unmodulated.tailIrradianceEnvelopeWm2, ...
    'candidateABER',aber,'candidateOP',outage,'candidateMeanPowerW',received, ...
    'statisticsMethod','Common fixed receiver frame; exact Rice tail; full realized phases; optimistic remote U-RIS tail envelope');
save(fullfile(folder,'comparison.mat'),'data','-v7.3');
end

function ids=choosePhase(primary,secondary,power,allowed)
% Exact OP=1 saturation is a tie, not evidence that a poor BER phase is best.
ids=zeros(1,size(primary,2));
for j=1:numel(ids)
    pool=find(allowed); best=min(primary(pool,j));
    pool=pool(primary(pool,j)<=best+8*eps(max(1,abs(best))));
    best=min(secondary(pool,j)); pool=pool(secondary(pool,j)<=best+8*eps(max(1,abs(best))));
    [~,k]=max(power(pool,j)); ids(j)=pool(k);
end
end
