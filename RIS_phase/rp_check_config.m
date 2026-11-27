function rp_check_config(c)
%RP_CHECK_CONFIG Check coupled physical and sampling parameters before simulation.
arguments
    c (1,1) struct
end
validateattributes(c.conditions,{'numeric'},{'real','finite','nonempty','ncols',2});
validateattributes(c.snapshotTimes,{'numeric'},{'real','finite','nonempty','vector','nonnegative'});
assert(all(c.snapshotTimes<=c.ris.duration),'rp:Times','Snapshot exceeds ris.duration.');
validateattributes([c.ris.Lx,c.ris.Ly,c.ris.thickness],{'numeric'},{'real','finite','positive'});
validateattributes([c.optics.nWater,c.optics.nGlass,c.optics.nAir],{'numeric'}, ...
    {'real','finite','positive'});
assert(c.optics.nGlass>c.optics.nWater && c.optics.nWater>=c.optics.nAir, ...
    'rp:Indices','Use nGlass > nWater >= nAir for the prescribed wetted glass.');
assert(c.optics.screenFraction>0 && c.optics.screenFraction<1,'rp:Screen','LC must be inside glass.');
validateattributes(c.array.count,{'numeric'},{'real','finite','size',[1 2],'integer','positive'});
validateattributes([c.array.pitchM;c.array.elementSizeM;c.array.apertureM], ...
    {'numeric'},{'real','finite','size',[3 2],'positive'});
assert(all(c.array.apertureM<=c.array.elementSizeM) && ...
    all(c.array.elementSizeM<=c.array.pitchM),'rp:Tiles','Require aperture <= tile size <= pitch.');
assert(all((c.array.count-1).*c.array.pitchM+c.array.elementSizeM<= ...
    [c.ris.Lx,c.ris.Ly]+1e-12),'rp:Tiles','Tile array does not fit the RIS board.');
validateattributes(c.source.wavelength,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(c.source.power,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(c.source.divergence,{'numeric'},{'scalar','real','finite','>',0,'<',pi/2});
validateattributes(c.source.edgeRelativeIntensity,{'numeric'},{'scalar','real','finite','>',0,'<',1});
validateattributes(c.device.maxAngleDeg,{'numeric'},{'scalar','real','finite','>',0,'<',90});
validateattributes(c.device.angleStepDeg,{'numeric'},{'scalar','real','finite','positive'});
assert(c.device.angleStepDeg<=c.device.maxAngleDeg,'rp:Step','Command step exceeds angle budget.');
assert(ismember(c.device.limitShape,{'box','disk'}),'rp:Limit','Choose box or disk.');
assert(strcmp(c.device.limitShape,'box'),'rp:Limit','Current device uses independent axes.');
assert(ismember(c.phase.model,{'linear','tanh','quadratic','quadratic_raw'}),'rp:PhaseModel','Unknown phase model.');
assert(~strcmp(c.phase.model,'tanh') || c.phase.separable,'rp:Integrability','Tanh requires separable=true.');
validateattributes(c.input.apertureOrder,{'numeric'},{'scalar','integer','>=',2});
validateattributes(c.mc.photons,{'numeric'},{'scalar','real','finite','integer','positive'});
validateattributes(c.mc.mapBins,{'numeric'},{'real','finite','size',[1 2],'integer','>=',8});
validateattributes(c.spot.mapBins,{'numeric'},{'real','finite','size',[1 2],'integer','>=',8});
validateattributes(c.spot.minIrradianceWm2,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(c.spot.coverageFraction,{'numeric'},{'scalar','real','finite','>',0,'<=',1});
validateattributes([c.device.electrodeWidthM,c.device.electrodeGapM,c.device.electrodePitchM], ...
    {'numeric'},{'real','finite','positive'});
assert(abs(c.device.electrodePitchM-c.device.electrodeWidthM-c.device.electrodeGapM)<1e-15, ...
    'rp:Electrodes','Electrode pitch must equal width plus gap.');
assert(abs(sind(c.device.maxAngleDeg)-c.source.wavelength/(2*c.device.electrodePitchM))<1e-12, ...
    'rp:ElectrodeLimit','Recompute the air angle reference from wavelength/(2*electrodePitchM).');
validateattributes(c.wave.windowPowerFraction,{'numeric'},{'scalar','real','finite','>',0,'<',1});
validateattributes(c.wave.windowPadding,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(c.wave.windowSamples,{'numeric'},{'scalar','integer','>=',3});
validateattributes(c.wave.apertureSamples,{'numeric'},{'scalar','real','finite','integer','>=',8});
validateattributes(c.wave.apertureOversampling,{'numeric'},{'scalar','real','finite','integer','>=',1});
validateattributes(c.wave.maximumOmittedPowerFraction,{'numeric'}, ...
    {'scalar','real','finite','>=',0,'<',1});
validateattributes(c.wave.fftPadding,{'numeric'},{'scalar','real','finite','integer','>=',2});
validateattributes(c.source.waist,{'numeric'},{'scalar','real','finite','positive'});
beamProduct=c.source.wavelength/(pi*c.optics.nWater);
assert(abs(c.source.waist*c.source.divergence-beamProduct)<=1e-8*beamProduct, ...
    'rp:GaussianBeam','The coherent Gaussian source requires waist*divergence=lambda0/(pi*nWater).');
validateattributes(c.detector.minimumRadiusPixels,{'numeric'},{'scalar','real','finite','>=',2});
validateattributes(c.detector.maxRefinements,{'numeric'},{'scalar','real','finite','integer','nonnegative'});
validateattributes(c.detector.zoomFactor,{'numeric'},{'scalar','real','finite','>',0,'<',1});
validateattributes(c.wave.maximumFrequencySamples,{'numeric'},{'scalar','real','finite','integer','>=',16});
validateattributes(c.wave.receiverSamplingMargin,{'numeric'},{'scalar','real','finite','>=',1});
validateattributes(c.detector.maximumUnresolvedPowerFraction,{'numeric'},{'scalar','real','finite','>=',0,'<',1});
validateattributes(c.rx.height,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(c.rx.apertureRadius,{'numeric'},{'scalar','real','finite','positive'});
end
