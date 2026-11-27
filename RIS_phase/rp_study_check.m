function rp_study_check(c)
%RP_STUDY_CHECK Validate units, receiver probabilities and finite-width bounds.
arguments
    c (1,1) struct
end
rp_check_config(c);
q=c.study;
validateattributes(q.distancesM,{'double'},{'size',[1 3],'positive','finite','increasing'});
validateattributes(q.widthsM,{'double'},{'vector','positive','finite','increasing'});
validateattributes(q.minimumWidthM,{'double'},{'scalar','positive','finite'});
validateattributes(q.maximumAxisRatio,{'double'},{'scalar','>=',1,'finite'});
validateattributes(q.directionSlackDeg,{'double'},{'scalar','nonnegative','finite','<',90});
validateattributes(q.positionSigmaM,{'double'},{'scalar','nonnegative','finite'});
validateattributes(q.angleSigmaRad,{'double'},{'scalar','nonnegative','finite'});
validateattributes(q.angleBiasRad,{'double'},{'size',[1 2],'real','finite'});
validateattributes(q.fovHalfAngleRad,{'double'},{'scalar','positive','finite','<',pi/2});
validateattributes([q.designPowerW,q.responsivityAW,q.noiseStdA,q.snrThreshold], ...
    {'double'},{'real','positive','finite'});
validateattributes(q.samples,{'double'},{'scalar','integer','>=',16});
validateattributes(q.validationSamples,{'double'},{'scalar','integer','>=',q.samples});
validateattributes(q.powerDBm,{'double'},{'vector','real','finite','increasing'});
validateattributes(q.positionSweepM,{'double'},{'vector','nonnegative','finite','increasing'});
validateattributes(q.angleSweepRad,{'double'},{'vector','nonnegative','finite','increasing'});
validateattributes(q.fovSweepRad,{'double'},{'vector','positive','finite','increasing','<',pi/2});
assert(all(q.widthsM>=q.minimumWidthM),'rp:WidthGrid','Requested widths must meet minimumWidthM.');
assert(ismember(q.turbulence,{'weak','strong','none'}),'rp:Turbulence','Use weak, strong or none.');
assert(ismember(q.atmosphere,{'weak','strong','none'}),'rp:Atmosphere','Use weak, strong or none GG.');
validateattributes([q.ggWeak,q.ggStrong],{'double'},{'size',[1 4],'positive','finite'});
validateattributes(q.robustExponent,{'double'},{'scalar','positive','finite'});
validateattributes(q.jitterGridM,{'double'},{'vector','nonnegative','finite','increasing'});
validateattributes(q.angleGridDeg,{'double'},{'vector','>=',0,'<',90,'finite','increasing'});
assert(isequal(q.families,{'gaussian','flat_top','robust'}),'rp:Families','Stage-two family set mismatch.');
assert(c.source.power==1,'rp:ReferencePower','The study uses a 1 W reference field; vary study.powerDBm instead.');
assert(isequal(c.array.count,[40 40]) && all(abs(c.array.pitchM-.025)<1e-12), ...
    'rp:DeviceGeometry','README requires 40x40 units at 25 mm center pitch.');
end
