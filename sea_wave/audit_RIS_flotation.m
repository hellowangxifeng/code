function audit = audit_RIS_flotation(p)
%AUDIT_RIS_FLOTATION Mass, displacement, both contact branches and dynamics.
% Accepts the actual simulation parameters; does not assume bottom z=0.
    if nargin==0, p=get_default_RIS_params(); end
    h=RIS_flotation(p);
    masses=[p.Lx*p.Ly*p.thickness*p.glassDensity, ...
        (prod([p.Lx p.Ly]+2*p.frameWidth)-p.Lx*p.Ly)*p.frameHeight*p.frameDensity, ...
        p.payloadMass,p.ballastMass];
    assert(abs(h.totalMassKg-sum(masses))<1e-10);
    % Independent root solve of the displaced-volume equation.
    volume=@(d) h.frameWaterplaneAreaM2*d+h.glassAreaM2*max(0,min(p.thickness,d-p.glassBottomAboveFrameBottom));
    draft=fzero(@(d) p.rho*volume(d)-sum(masses),[0 p.frameHeight]);
    assert(abs(h.frameDraftM-draft)<1e-11);
    assert(abs(h.buoyancyResidualN)<1e-10 && all(h.restoring>0));
    assert(abs(h.glassTopZ-h.glassBottomZ-p.thickness)<1e-12);
    assert(all(h.inertiaKgM2>0));

    % Independently test the dry/wet slope change at exact bottom contact.
    % Move only the TEST fixture's mounting height to its contact waterline;
    % the production 1 m plate and its fixed mounting height remain unchanged.
    contact=p;
    contact.glassBottomAboveFrameBottom=h.totalMassKg/(p.rho*h.frameWaterplaneAreaM2);
    assert(contact.glassBottomAboveFrameBottom+p.thickness<p.frameHeight);
    hc=RIS_flotation(contact);
    assert(abs(hc.glassBottomZ)<1e-12);
    lighter=contact; lighter.ballastMass=contact.ballastMass-0.1;
    heavier=contact; heavier.ballastMass=contact.ballastMass+0.1;
    hl=RIS_flotation(lighter); hh=RIS_flotation(heavier);
    assert(abs(hl.glassBottomZ-0.1/(p.rho*h.frameWaterplaneAreaM2))<1e-12);
    assert(abs(hh.glassBottomZ+0.1/(p.rho*prod(h.outerSizeM)))<1e-12);

    frequency=0.5; omega=2*pi*frequency; k=omega^2/p.g;
    cfg=p; cfg.duration=0.2; cfg.initialization='rest'; cfg.radiationFrequencyBins=16;
    cfg.initialPosition=[]; cfg.initialEulerAngles=[]; cfg.initialBodyVelocity=[];
    cfg.horizontalInitialPosition=[0;0]; cfg.horizontalInitialVelocity=[];
    cfg.wave=struct('omega',omega,'kx',0.6*k,'ky',0.8*k,'coefficient',0);
    calm=RIS_float(cfg);
    assert(max(abs(calm.eta),[],'all')<1e-12);
    assert(max(abs(calm.bottomCenter(3,:)-h.glassBottomZ))<1e-12);
    assert(max(abs(calm.boardCenter(3,:)-h.glassTopZ))<1e-12);
    assert(max(abs(calm.rigidBodyMassMatrix(1:3,1:3)-h.totalMassKg*eye(3)),[],'all')<1e-12);
    cfg.wave.coefficient=1e-4;
    regular=RIS_float(cfg);
    assert(regular.validity.numericalSolutionFinite);
    assert(regular.validity.maxEquationResidual<1e-9);
    assert(regular.validity.maxRotationOrthogonalityError<1e-12);
    separation=regular.boardCenter-regular.bottomCenter-p.thickness*regular.normal;
    assert(max(abs(separation),[],'all')<1e-12);
    audit=struct('passed',true,'glassSizeM',[p.Lx p.Ly], ...
        'flotation',h,'independentDraftM',draft,'contactTestMountM',contact.glassBottomAboveFrameBottom, ...
        'massPerturbationKg',0.1,'lighterBottomZM',hl.glassBottomZ,'heavierBottomZM',hh.glassBottomZ, ...
        'calmMaxMotion',max(abs(calm.eta),[],'all'), ...
        'regularWaveEquationResidual',regular.validity.maxEquationResidual);
end
