function d=rp_stage2_select(c,d)
%RP_STAGE2_SELECT Compare actual maps under the stage-two channel model.
rows=zeros(numel(d.bank),11);
for j=1:numel(d.bank)
    b=d.bank{j}; f=b.spot; m=b.designMetrics;
    tilt=atan2d(norm(f.centerXY-d.geometry.pRis(1:2).'),d.distanceM);
    rows(j,:)=[j,b.targetWidthM,100*f.equivalentRadiusM,f.minorRadiusM, ...
        f.axisRatio,tilt,m.meanReceivedPowerW,m.aber,m.op,f.r50M,f.r90M];
end
T=array2table(rows,'VariableNames',{'candidate','targetWidthM','actualRadiusCM', ...
    'minorRadiusM','axisRatio','tiltDeg','meanReceivedPowerW','aber','op','r50M','r90M'});
T.family=string(cellfun(@(b) b.family,d.bank,'UniformOutput',false)).';
feasible=T.minorRadiusM>=c.study.minimumWidthM & T.actualRadiusCM>=100*c.study.minimumWidthM & ...
    T.axisRatio<=c.study.maximumAxisRatio & T.r50M>=.02 & ismember(T.family,["robust","flat_top"]);
assert(any(feasible),'rp:Stage2Feasible','No finite robust/flat-top design.');
feasible=feasible & T.tiltDeg<=min(T.tiltDeg(feasible))+c.study.directionSlackDeg;
T.feasible=feasible;
scores=[T.aber,T.op,-T.meanReceivedPowerW]; scores(~feasible,:)=Inf;
% Keep the robust and flat-top families independent, so the adaptive flat-top
% comparison can genuinely win in a sigma interval rather than being a subset.
robust=feasible & T.family=="robust";
assert(any(robust),'rp:RobustFeasible','No finite robust-family design.');
scores(~robust,1:2)=Inf;
[~,selected]=min(scores,[],1); d.objectiveIndices=selected;
ids=find(T.family=="flat_top"); [~,nearest]=min(abs(T.targetWidthM(ids)-c.study.referenceWidthM));
flat=ids(nearest);
ids=find(T.family=="gaussian"); [~,nearest]=min(abs(T.targetWidthM(ids)-c.study.referenceWidthM));
d.gaussianIndex=ids(nearest);
d.comparisonIndices=[1 2 flat selected];
d.comparisonNames={'Unmodulated RIS','Steering-only','Flat top','Robust ABER','Robust OP','Max mean power'};
d.objectiveNames={'ABER','OP','Power'}; d.summary=T;
end
