function result=rp_stage2_angle_case(c,d,angleId)
%RP_STAGE2_ANGLE_CASE Independent resumable incidence case for dynamic scheduling.
rp_setup(); angle=c.study.angleGridDeg(angleId);
folder=fullfile(c.output,'incidence_cases'); if ~isfolder(folder), mkdir(folder); end
file=fullfile(folder,sprintf('z%g_a%05.2f.mat',d.distanceM,angle));
spec=cellfun(@(b) struct('family',b.family,'width',b.targetWidthM,'margin',b.margin), ...
    d.bank(d.comparisonIndices(3:6)));
model=[c.study.signature,jsonencode(spec),fileread([mfilename('fullpath') '.m']), ...
    fileread(which('rp_direct_wave')),fileread(which('rp_direct_channel')),fileread(which('rp_direct_density'))];
if isfile(file)
    cached=load(file,'result','model');
    if strcmp(model,cached.model), result=cached.result; return; end
end
[g,mc,a,cc]=rp_stage2_scene(c,d.distanceM,angle);
bank=cell(1,6); bank{1}=rp_direct_wave(cc,g);
bank{2}=rp_stage2_wave(cc,g,mc,a,'steering',0,0);
for k=3:6
    nominal=d.bank{d.comparisonIndices(k)};
    bank{k}=rp_stage2_wave(cc,g,mc,a,nominal.family,nominal.targetWidthM,nominal.margin);
end
rows=zeros(6,8);
for k=1:6
    b=bank{k}; m=b.designMetrics;
    rows(k,:)=[angleId,angle,k,m.meanReceivedPowerW,m.aber,m.op,b.spot.equivalentRadiusM,b.spot.axisRatio];
end
result=struct('table',array2table(rows,'VariableNames',{'angleId','incidenceDeg','schemeId', ...
    'meanReceivedPowerW','aber','op','equivalentRadiusM','axisRatio'}),'bank',{bank}, ...
    'distanceM',d.distanceM,'angleId',angleId,'incidenceDeg',angle);
save(file,'result','model','-v7.3');
fprintf('Incidence z=%g angle=%g complete\n',d.distanceM,angle);
end
