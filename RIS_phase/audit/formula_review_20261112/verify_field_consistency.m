function review=verify_field_consistency(s)
%VERIFY_FIELD_CONSISTENCY Reassemble every stored scatter map from raw packets.
arguments
    s (1,1) struct
end
c=s.config;
assert(strcmp(rp_study_signature(c),c.study.signature),'rp:ReviewSignature', ...
    'Numerical sources differ from the saved study. Reassess propagation before reusing figures.');
rows=zeros(sum(cellfun(@(d) numel(d.bank),s.distances)),7);
geometryRows=zeros(numel(s.distances),5); row=0;
for j=1:numel(s.distances)
    d=s.distances{j}; g=d.geometry; a=d.illumination;
    expectedSource=g.RRis.'*(g.pTx-g.pRis);
    sourceError=max(abs(expectedSource-g.sourceLocal));
    rebuilt=rp_illuminate(c,g,d.mc);
    lcError=max(abs(rebuilt.lcXY-a.lcXY),[],'all');
    powerError=max(abs(rebuilt.powerW-a.powerW));
    assert(sourceError<1e-13 && lcError<1e-13 && powerError<1e-13,'rp:ReviewLC', ...
        'Saved local source or bottom-to-LC transport differs from current implementation.');
    geometryRows(j,:)=[d.distanceM,g.pRis(3),sourceError,lcError,powerError];
    for k=1:numel(d.bank)
        item=d.bank{k}; w=item.wave;
        ray=rp_trace(c,g,a,item.phase);
        use=ray.valid & ~a.ballistic;
        dx=w.x(2)-w.x(1); dy=w.y(2)-w.y(1);
        xe=[w.x-dx/2,w.x(end)+dx/2];
        ye=[w.y-dy/2,w.y(end)+dy/2];
        ix=discretize(ray.landing(use,1),xe);
        iy=discretize(ray.landing(use,2),ye);
        power=ray.powerW(use); inside=isfinite(ix) & isfinite(iy);
        halo=accumarray([iy(inside),ix(inside)],power(inside), ...
            [numel(w.y),numel(w.x)],@sum,0)/w.cellAreaM2;
        haloError=max(abs(halo-w.scatteredIrradianceWm2),[],'all');
        sumError=max(abs(w.irradianceWm2-w.ballisticPixelMeanWm2-halo),[],'all');
        areaError=abs(w.scatteredMapPowerW-sum(power(inside)));
        directWeightScale=max(abs(a.ballisticModelScale-1));
        assert(haloError<1e-11 && sumError<1e-11 && areaError<1e-12, ...
            'rp:ReviewHalo','Saved scatter intensity or hybrid sum disagrees with packet transport.');
        assert(directWeightScale==0,'rp:ReviewScale','Unexpected fitted ballistic normalization.');
        row=row+1;
        rows(row,:)=[d.distanceM,k,nnz(use),haloError,sumError,areaError,w.scatteredMapPowerW];
    end
end
review=struct('status','passed','numericalSignature',c.study.signature, ...
    'candidateCount',row,'geometry',array2table(geometryRows,'VariableNames', ...
    {'distanceM','boardWorldZM','sourceLocalErrorM','lcPositionErrorM','packetPowerErrorW'}), ...
    'scattering',array2table(rows,'VariableNames', ...
    {'distanceM','candidate','validScatterPackets','haloMaxErrorWm2','sumMaxErrorWm2', ...
    'scatterPowerErrorW','scatterMapPowerW'}));
folder=fileparts(mfilename('fullpath'));
writetable(review.geometry,fullfile(folder,'geometry_checks.csv'));
writetable(review.scattering,fullfile(folder,'all_candidate_scatter_checks.csv'));
save(fullfile(folder,'field_consistency.mat'),'review');
fprintf('Field consistency passed for %d stored candidates at %d distances.\n',row,numel(s.distances));
end
