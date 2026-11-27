function result=rp_stage2_merge_angles(c,cases,z)
%RP_STAGE2_MERGE_ANGLES Gather disjoint incidence cases without concurrent writers.
maps=cell(numel(c.study.angleGridDeg),6); T=[];
for j=1:numel(cases)
    q=cases{j}; if abs(q.distanceM-z)>1e-9, continue; end
    T=[T;q.table]; maps(q.angleId,:)=q.bank; %#ok<AGROW>
end
T=sortrows(T,{'angleId','schemeId'});
assert(height(T)==6*numel(c.study.angleGridDeg) && all(~cellfun(@isempty,maps),'all'), ...
    'rp:IncidenceCases','Missing or duplicate incidence cases.');
result=struct('table',T,'maps',{maps},'distanceM',z, ...
    'criticalAngleDeg',asind(c.optics.nAir/c.optics.nWater));
save(fullfile(c.output,sprintf('incidence_z%g.mat',z)),'result','-v7.3');
writetable(T,fullfile(c.output,sprintf('incidence_z%g.csv',z)));
end
