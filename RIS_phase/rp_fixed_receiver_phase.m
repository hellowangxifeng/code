function p=rp_fixed_receiver_phase(c,d,family,width,margin)
%RP_FIXED_RECEIVER_PHASE One geometric receiver target, independent of metrics.
g=d.geometry; a=d.illumination;
p=rp_finite_phase(c,g,a,'flat_top',.10,margin);
target=p.targetCenterXY;
if ismember(family,{'flat_top','robust'})
    p=rp_finite_phase(c,g,a,family,width,margin);
elseif strcmp(family,'steering')
    % One constant LC ramp. No curvature compensation or focusing pistons.
    % Aiming the ramp at the common RX keeps steering-only a fair baseline.
    p=d.bank{2}.phase; weights=a.unitDirectPowerW(p.unitIds); weights=weights/sum(weights);
    reference=sum(p.referenceXY.*weights,1);
    incoming=rp_incident_field(c,g,reference);
    outgoing=rp_paraxial(c,g,reference,target);
    [b,command]=rp_command(c,outgoing.q-incoming.q);
    p.gradientOverK0=repmat(b,numel(p.unitIds),1);
    p.commandDeg=repmat(command,numel(p.unitIds),1); p.phaseCurvature(:)=0;
    p.phaseAtReferenceRad=mod(2*pi/c.source.wavelength*(p.referenceXY*b.'),2*pi);
    p.family='steering'; p.phaseModel='linear'; p.targetWidthM=0;
    actual=zeros(numel(p.unitIds),2); valid=false(numel(p.unitIds),1);
    for k=1:numel(p.unitIds)
        [actual(k,:),valid(k)]=rp_carrier_landing(c,g,p.referenceXY(k,:),p.qInAtReference(k,:)+b);
    end
    p.targetsXY=actual; p.outgoing=rp_paraxial(c,g,p.referenceXY,actual);
    p.outgoing.valid=p.outgoing.valid & valid;
    p.objective='constant global phase ramp toward a fixed receiver; no reshaping';
else
    error('rp:FixedFamily','Unsupported fixed-target design family.');
end
p.targetCenterXY=target; p.targetMargin=margin;
end
