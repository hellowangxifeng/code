function p=rp_zero_phase(c,g,a)
%RP_ZERO_PHASE Unmodulated physical RIS: same glass, masks, zero LC phase.
p=rp_phase(c,g,a); p.gradientOverK0(:)=0; p.commandDeg(:)=0;
p.phaseCurvature(:)=0; p.phaseAtReferenceRad(:)=0;
p.phaseModel='linear'; p.family='unmodulated'; p.targetWidthM=0;
q=p.qInAtReference; valid=sum(q.^2,2)<c.optics.nAir^2;
targets=repmat(g.pRis(1:2).',numel(p.unitIds),1);
for j=find(valid).'
    [targets(j,:),valid(j)]=rp_carrier_landing(c,g,p.referenceXY(j,:),q(j,:));
end
p.targetsXY=targets; p.outgoing=rp_paraxial(c,g,p.referenceXY,targets);
p.outgoing.valid=p.outgoing.valid & valid;
p.referenceDirectionWorld=[q,sqrt(max(0,c.optics.nAir^2-sum(q.^2,2)))]/c.optics.nAir*g.RRis.';
p.objective='zero physical LC phase, identical tiled transparent apertures';
end
