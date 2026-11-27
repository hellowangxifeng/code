function out = RIS_float(cfg)
%RIS_FLOAT Glass RIS + load-bearing frame, six-DOF reference hydrodynamics.
% All angles now use right-hand ZYX Euler conventions (old pitch sign changed).
% eta=[surge;sway;heave;roll;pitch;yaw], bodyVelocity=[u;v;w;p;q;r].
% q=[heave;pitch;roll;surge;sway;yaw] retains the first five ROW positions.
% Exact quaternion kinematics, rigid gyroscopic and reference-added-inertia coupling.
% Radiation memory is a passive fit to the ORIGINAL vertical BEM, with errors.
% Incident forces are linearized about the initial horizontal heading and
% reference position. Horizontal/yaw diffraction, radiation and viscosity are
% NOT solved. This is a six-DOF approximate hydrodynamic model, not validation
% of complete six-DOF hydrodynamic coefficients. See RIS_float_model.md.
% cfg.wave accepts original oceanSurfaceModel components without reseeding.
% This function is compute-only; figures are owned by run_RIS_wind_comparison.
% cfg.initialization='rest' starts all states and radiation memory at rest.
% initialPosition/EulerAngles/bodyVelocity specify initial conditions only.
    arguments
        cfg (1,1) struct = struct()
    end
    p = get_default_RIS_params();
    fields = fieldnames(cfg);
    for j = 1:numel(fields)
        if ~isfield(p,fields{j})
            error('RIS_float:Option','Unknown option: %s. See get_default_RIS_params.',fields{j});
        end
        p.(fields{j}) = cfg.(fields{j});
    end
    p = validateParameters(p);
    flotation = RIS_flotation(p);
    draft = flotation.frameDraftM;
    zG = flotation.centerOfMassZ;
    rigidMass = [repmat(flotation.totalMassKg,3,1);flotation.inertiaKgM2];
    restoring = flotation.restoring;
    % Thin-waterline BEM on the load-bearing ring, not on the optical glass.
    % A statically submerged glass activates the filled reference footprint.
    hydro = p;
    hydro.Lx = flotation.outerSizeM(1); hydro.Ly = flotation.outerSizeM(2);
    hydro.innerLx = p.Lx; hydro.innerLy = p.Ly;
    if flotation.glassDraftM>1e-12
        hydro.innerLx = 0; hydro.innerLy = 0;
    end
    spacing = [hydro.Lx,hydro.Ly]/p.panelCount;
    % Only a ring reference has inner boundaries to align. A statically wet
    % glass uses the filled rectangle, so imposing the old ring test is wrong.
    if hydro.innerLx>0 && ...
            any(abs([p.frameWidth,p.frameWidth]./spacing-round([p.frameWidth,p.frameWidth]./spacing))>1e-8)
        error('RIS_float:FrameMesh','panelCount must align panel edges with the inner frame edges.');
    end
    [wave,spectrum] = makeWave(p);
    heading = 0;
    if ~isempty(p.initialEulerAngles), heading = p.initialEulerAngles(3); end
    [headingR,~] = RIS_rotation([0;0;heading]);
    bodyWave = wave;
    bodyK = headingR(1:2,1:2).'*[wave.kx.';wave.ky.'];
    bodyWave.kx = bodyK(1,:).'; bodyWave.ky = bodyK(2,:).';
    fitFrequency = 2*pi*logspace(log10(p.fMinHz),log10(p.fMaxHz),p.radiationFrequencyBins).';
    frequency = unique([wave.omega;fitFrequency]);
    [~,index] = ismember(wave.omega,frequency);
    % BEM modes are [heave, x-height slope, y-height slope]; pitch=-x-slope.
    modalMass = rigidMass([3,5,4]);
    modalK = restoring([1,3,2]);
    [oldRao,A,B,F] = plateResponse(hydro,frequency,bodyWave.kx,bodyWave.ky,index,modalMass,modalK);
    radiation = RIS_radiation(frequency,A,B,p.radiationModes);
    % Convert all BEM results to [heave, roll, pitch].
    modes = [3,4,5];
    radiation.residues = radiation.residues(:,[1,3,2]);
    radiation.infiniteAddedMass = radiation.infiniteAddedMass([1,3,2]);
    radiation.relativeImpedanceRmsError = radiation.relativeImpedanceRmsError([1,3,2]);
    radiation.fittedAddedMass = radiation.fittedAddedMass([1,3,2],:);
    radiation.fittedDamping = radiation.fittedDamping([1,3,2],:);
    totalMass = rigidMass;
    totalMass(modes) = totalMass(modes)+radiation.infiniteAddedMass;
    % Outer walls minus open inner walls. The submerged glass closes the
    % upper part of the inner opening; these pressure integrals retain yaw.
    excitation = RIS_sidewall_load(hydro,bodyWave,draft,zG)-RIS_sidewall_load(p,bodyWave,draft,zG);
    if flotation.glassDraftM>1e-12
        excitation = excitation+RIS_sidewall_load(p,bodyWave,flotation.glassDraftM,zG);
    end
    excitation(3,:) = F(1,:);
    excitation(4,:) = excitation(4,:)+F(3,:);
    excitation(5,:) = excitation(5,:)-F(2,:);
    initialXY = p.horizontalInitialPosition;
    if ~isempty(p.initialPosition), initialXY = p.initialPosition(1:2); end
    referenceCenter = p.center+initialXY;
    c = wave.coefficient.*exp(1i*(wave.omega*p.startTime- ...
        wave.kx*referenceCenter(1)-wave.ky*referenceCenter(2)));
    forceHat = excitation.*c.';
    model = struct('mass',totalMass,'rigidMass',rigidMass,'restoring',restoring, ...
        'radiation',radiation,'headingR',headingR,'forceHat',forceHat,'omega',wave.omega);
    [state,admittance] = initialState(p,model,initialXY,heading);
    rao = admittance.*excitation;
    t = unique([0:p.dt:p.duration,p.duration]);
    [states,step] = integrate(model,state,t,p.integrationStep);
    out = assembleOutput(p,model,t,states,wave,spectrum,referenceCenter,flotation,step);
    out.freeboard = flotation.glassTopZ;
    out.flotation = flotation;
    out.addedMass = A;
    out.radiationDamping = B;
    out.hydroFrequencyHz = frequency/(2*pi);
    out.excitationPerMetre = excitation;
    out.responseAmplitudeOperator = [oldRao(1,:);-oldRao(2,:);oldRao(3,:)];
    out.responseAmplitudeOperator6 = rao;
    out.responseComponents = rao.*c.';
end

function [state,rao] = initialState(p,m,initialXY,heading)
    s = 1i*m.omega.';
    rad = m.radiation;
    basis = (s+rad.width)./((s+rad.width).^2+rad.frequency.^2);
    impedance = rad.residues.'*basis;
    stiffness = [0;0;m.restoring;0];
    Z = m.mass.*s.^2+stiffness;
    Z([3,4,5],:) = Z([3,4,5],:)+s.*impedance;
    response = m.forceHat./Z;
    rao = 1./Z;
    state = zeros(13+6*numel(rad.width),1);
    state(1:3) = [initialXY;0];
    angles = [0;0;heading];
    if strcmp(p.initialization,'waveSteady')
        displacement = real(sum(response,2));
        velocity = real(sum(response.*s,2));
        state(3) = displacement(3);
        angles(1:2) = displacement(4:5);
        worldVelocity = m.headingR*velocity(1:3);
        angular = velocity(4:6); angular(3) = 0;
        [rotation,~] = RIS_rotation(angles);
        state(8:10) = rotation.'*worldVelocity;
        state(11:13) = angular;
        modalVelocity = response([3,4,5],:).*s;
        memoryX = real(basis*modalVelocity.');
        memoryY = real((rad.frequency./((s+rad.width).^2+rad.frequency.^2))*modalVelocity.');
        state(14:end) = [memoryX(:);memoryY(:)];
    end
    if ~isempty(p.initialPosition), state(1:3) = p.initialPosition; end
    if ~isempty(p.initialEulerAngles), angles = p.initialEulerAngles; end
    [rotation,quaternion] = RIS_rotation(angles);
    state(4:7) = quaternion;
    if ~isempty(p.horizontalInitialVelocity)
        worldVelocity = rotation*state(8:10);
        worldVelocity(1:2) = p.horizontalInitialVelocity;
        state(8:10) = rotation.'*worldVelocity;
    end
    if ~isempty(p.initialBodyVelocity), state(8:13) = p.initialBodyVelocity; end
end

function [states,largestStep] = integrate(m,state,t,requestedStep)
% Fixed-step RK4 with stability subdivisions; dt only controls saved samples.
    poleRate = max(hypot(m.radiation.width,m.radiation.frequency));
    naturalRate = sqrt(max(m.restoring./m.mass([3,4,5])));
    safeStep = min([requestedStep,0.6/poleRate,0.3/naturalRate]);
    states = zeros(numel(state),numel(t)); states(:,1) = state;
    largestStep = 0;
    for j = 2:numel(t)
        steps = ceil((t(j)-t(j-1))/safeStep);
        h = (t(j)-t(j-1))/steps;
        largestStep = max(largestStep,h);
        % Harmonic excitation evaluated exactly at every RK stage.
        stageTimes = t(j-1)+(0:2*steps)*h/2;
        forces = real(m.forceHat*exp(1i*(m.omega*stageTimes)));
        for k = 1:steps
            first = 2*k-1;
            k1 = derivative(state,m,forces(:,first));
            k2 = derivative(state+h/2*k1,m,forces(:,first+1));
            k3 = derivative(state+h/2*k2,m,forces(:,first+1));
            k4 = derivative(state+h*k3,m,forces(:,first+2));
            state = state+h/6*(k1+2*k2+2*k3+k4);
            state(4:7) = state(4:7)/norm(state(4:7));
        end
        if any(~isfinite(state))
            error('RIS_float:Integration','Nonfinite state at t=%g; reduce integrationStep.',t(j));
        end
        states(:,j) = state;
    end
end

function [d,load] = derivative(state,m,force)
% Exact rigid-body equations in world axes; linear hydrodynamics in REF axes.
    [R,q] = RIS_rotation(state(4:7));
    nu = state(8:13); v = nu(1:3); omega = nu(4:6);
    ref = m.headingR;
    worldVelocity = R*v; worldOmega = R*omega;
    worldMass = diag([m.rigidMass(1),m.rigidMass(2),m.mass(3)]);
    worldInertia = R*diag(m.rigidMass(4:6))*R.';
    addedInertia = ref*diag([m.radiation.infiniteAddedMass(2:3);0])*ref.';
    effectiveInertia = worldInertia+addedInertia;
    effectiveBody = blkdiag(R.'*worldMass*R,R.'*effectiveInertia*R);
    coupling = [effectiveBody(1:3,1:3)*cross(omega,v);cross(omega,m.rigidMass(4:6).*omega)];
    n = R(3,:).';
    grad = [m.restoring(3)*n(1);m.restoring(2)*n(2);0];
    hydrostatic = [R.'*[0;0;-m.restoring(1)*state(3)];cross(n,grad)];
    waveLoad = [R.'*ref*force(1:3);R.'*ref*force(4:6)];
    nr = numel(m.radiation.width);
    x = reshape(state(14:13+3*nr),nr,3);
    y = reshape(state(14+3*nr:end),nr,3);
    memory = sum(m.radiation.residues.*x,1).';
    radiationLoad = [R.'*ref*[0;0;-memory(1)];R.'*ref*[-memory(2:3);0]];
    acceleration = effectiveBody\(waveLoad+hydrostatic+radiationLoad-coupling);
    refVelocity = [worldVelocity(3);ref(:,1:2).'*worldOmega];
    dx = -m.radiation.width.*x-m.radiation.frequency.*y+refVelocity.';
    dy = m.radiation.frequency.*x-m.radiation.width.*y;
    dq = 0.5*[-q(2:4).'*omega;q(1)*omega+cross(q(2:4),omega)];
    d = [worldVelocity;dq;acceleration;dx(:);dy(:)];
    if nargout>1
        load = struct('wave',waveLoad,'hydrostatic',hydrostatic,'radiation',radiationLoad, ...
            'coriolis',coupling,'acceleration',acceleration,'effectiveMass',effectiveBody, ...
            'residual',effectiveBody*acceleration+coupling-waveLoad-hydrostatic-radiationLoad);
    end
end

function out = assembleOutput(p,m,t,states,wave,spectrum,referenceCenter,h,step)
    count = numel(t);
    rotation = zeros(3,3,count); angles = zeros(3,count); vel = zeros(3,count);
    effectiveMass = zeros(6,6,count);
    acceleration = zeros(6,count); loads = zeros(6,count,4); residual = zeros(6,count);
    corners = zeros(4,3,count); center = states(1:3,:)+[p.center;h.centerOfMassZ];
    local = [-p.Lx/2,-p.Ly/2,h.glassTopFromCGM;p.Lx/2,-p.Ly/2,h.glassTopFromCGM; ...
        p.Lx/2,p.Ly/2,h.glassTopFromCGM;-p.Lx/2,p.Ly/2,h.glassTopFromCGM].';
    normal = zeros(3,count); boardCenter = normal; bottomCenter = normal; eulerRate = normal;
    forces = real(m.forceHat*exp(1i*(m.omega*t)));
    energy = zeros(1,count); orthogonality = zeros(1,count); power = zeros(1,count);
    for j = 1:count
        [R,~,angles(:,j)] = RIS_rotation(states(4:7,j));
        rotation(:,:,j) = R;
        normal(:,j) = R(:,3);
        boardCenter(:,j) = center(:,j)+R(:,3)*h.glassTopFromCGM;
        bottomCenter(:,j) = center(:,j)+R(:,3)*h.glassBottomFromCGM;
        corners(:,:,j) = (R*local+center(:,j)).';
        vel(:,j) = R*states(8:10,j);
        a = angles(1,j); b = angles(2,j);
        T = [1,sin(a)*tan(b),cos(a)*tan(b);0,cos(a),-sin(a); ...
            0,sin(a)/cos(b),cos(a)/cos(b)];
        eulerRate(:,j) = T*states(11:13,j);
        [~,load] = derivative(states(:,j),m,forces(:,j));
        acceleration(:,j) = load.acceleration;
        effectiveMass(:,:,j) = load.effectiveMass;
        loads(:,j,:) = reshape([load.wave,load.hydrostatic,load.radiation,load.coriolis],6,1,4);
        residual(:,j) = load.residual;
        nr = numel(m.radiation.width);
        x = reshape(states(14:13+3*nr,j),nr,3); y = reshape(states(14+3*nr:end,j),nr,3);
        memoryEnergy = 0.5*sum(m.radiation.residues.*(x.^2+y.^2),'all');
        potential = 0.5*(m.restoring(1)*states(3,j)^2+ ...
            m.restoring(2)*R(3,2)^2+m.restoring(3)*R(3,1)^2);
        energy(j) = 0.5*states(8:13,j).'*load.effectiveMass*states(8:13,j)+potential+memoryEnergy;
        power(j) = states(8:13,j).'*load.wave- ...
            sum(m.radiation.width.*m.radiation.residues.*(x.^2+y.^2),'all');
        orthogonality(j) = norm(R.'*R-eye(3),'fro');
    end
    angles(3,:) = unwrap(angles(3,:));
    eta = [states(1:3,:);angles]; order = [3,5,4,1,2,6];
    q = eta(order,:); qd = [vel;eulerRate]; qd = qd(order,:);
    field = struct('wave',wave,'parameters',p);
    etaCenter = RIS_wave_field(field,center(1:2,:).',t).';
    tilt = acosd(max(-1,min(1,normal(3,:))));
    k = wave.omega.^2/p.g;
    orbitalWeight = sum(abs(m.forceHat(1:2,:)./(m.rigidMass(1)*wave.omega.'.^2)).^2,1).';
    shortFraction = sum(orbitalWeight(k*max(p.Lx,p.Ly)>pi))/max(sum(orbitalWeight),realmin);
    % The original encounter-frequency diagnostic, without a waves-by-times array.
    phaseGradient = [wave.kx./wave.omega,wave.ky./wave.omega].*abs(wave.coefficient);
    encounterCovariance = phaseGradient.'*phaseGradient;
    dopplerVariance = sum(vel(1:2,:).*(encounterCovariance*vel(1:2,:)),1);
    encounterRatio = sqrt(max(0,mean(dopplerVariance))/max(sum(abs(wave.coefficient).^2),realmin));
    % Sample a 5x5 footprint; a rigid plane cannot match a curved sea surface.
    [sx,sy] = ndgrid(linspace(-p.Lx/2,p.Lx/2,5),linspace(-p.Ly/2,p.Ly/2,5));
    bottomGap = zeros(numel(sx),count); topClearance = bottomGap;
    for j = 1:numel(sx)
        point = bottomCenter+reshape(rotation(:,1,:),3,[])*sx(j)+ ...
            reshape(rotation(:,2,:),3,[])*sy(j);
        seaBottom = RIS_wave_field(field,point(1:2,:).',t).';
        topPoint = point+p.thickness*normal;
        seaTop = RIS_wave_field(field,topPoint(1:2,:).',t).';
        bottomGap(j,:) = point(3,:)-seaBottom;
        topClearance(j,:) = topPoint(3,:)-seaTop;
    end
    mismatchByTime = max(abs(bottomGap),[],1);
    mismatch = max(mismatchByTime);
    validity = struct('maxTiltDeg',max(tilt),'smallAngleUnder10Deg',all(tilt<=10), ...
        'maxHeadingChangeDeg',max(abs(rad2deg(angles(3,:)-angles(3,1)))), ...
        'hydrodynamicsComputedWithBEM',true,'usesExternalHydrodynamicDatabase',false, ...
        'verticalRadiationDiffractionSolved',true,'radiationMemoryUsesFit',true, ...
        'horizontalDiffractionSolved',false,'yawDiffractionSolved',false, ...
        'horizontalYawRadiationSolved',false,'fullSecondOrderDriftSolved',false, ...
        'directWindLoadIncluded',false,'experimentallyCalibrated',false, ...
        'horizontalShortWaveVarianceFraction',shortFraction, ...
        'rmsEncounterFrequencyChange',encounterRatio, ...
        'maxIncidentWaveBoardMismatchM',mismatch,'incidentMismatchExceedsDraft',mismatch>h.frameDraftM, ...
        'glassWetDryContactSolved',false,'frameThinWaterlineApproximation',true, ...
        'hydrodynamicReferenceIsRing',h.glassDraftM<=1e-12, ...
        'bottomGapSampleFraction',mean(bottomGap>1e-4,'all'), ...
        'topWettedSampleFraction',mean(topClearance<0,'all'), ...
        'minTopClearanceM',min(topClearance,[],'all'), ...
        'numericalSolutionFinite',all(isfinite(states(:))), ...
        'maxRotationOrthogonalityError',max(orthogonality), ...
        'maxEquationResidual',max(abs(residual),[],'all'), ...
        'radiationFitRelativeRmsError',m.radiation.relativeImpedanceRmsError, ...
        'referenceLinearHydrodynamics',true,'meshConvergenceCheckedThisCall',false);
    out = struct('t',t,'absoluteTime',p.startTime+t,'eta',eta,'q',q, ...
        'dofOrder',{{'surge','sway','heave','roll','pitch','yaw'}}, ...
        'qOrder',{{'heave','pitch','roll','surge','sway','yaw'}}, ...
        'surge',q(4,:),'sway',q(5,:),'heave',q(1,:), ...
        'roll',q(3,:),'pitch',q(2,:),'yaw',q(6,:), ...
        'rollDeg',rad2deg(q(3,:)),'pitchDeg',rad2deg(q(2,:)),'yawDeg',rad2deg(q(6,:)), ...
        'velocity',qd,'bodyVelocity',states(8:13,:),'bodyAcceleration',acceleration, ...
        'rotationMatrix',rotation,'quaternion',states(4:7,:),'normal',normal, ...
        'corners',corners,'boardCenter',boardCenter,'bottomCenter',bottomCenter,'centerOfMass',center, ...
        'bottomSeaGapM',bottomGap,'topSeaClearanceM',topClearance, ...
        'centerTrajectory',center(1:2,:),'horizontalVelocity',vel(1:2,:), ...
        'waveLoadBody',loads(:,:,1),'hydrostaticLoadBody',loads(:,:,2), ...
        'radiationLoadBody',loads(:,:,3),'coriolisLoadBody',loads(:,:,4), ...
        'waveLoadReference',forces,'horizontalWaveForce',m.headingR(1:2,1:2)*forces(1:2,:), ...
        'yawWaveMoment',forces(6,:),'equationResidual',residual, ...
        'energy',energy,'netPower',power,'integrationStepUsed',step, ...
        'mass',m.rigidMass([3,5,4]),'massMatrix',diag(m.rigidMass(order)), ...
        'bodyMassMatrix',effectiveMass,'rigidBodyMassMatrix',diag(m.rigidMass), ...
        'stiffness',[m.restoring([1,3,2]);0;0;0], ...
        'stiffnessMatrix',diag([m.restoring([1,3,2]);0;0;0]), ...
        'radiationModel',m.radiation,'etaCenter',etaCenter, ...
        'incidentWaveBoardMismatchM',mismatchByTime, ...
        'etaReferenceCenter',real((wave.coefficient.*exp(1i*(wave.omega*p.startTime- ...
            wave.kx*referenceCenter(1)-wave.ky*referenceCenter(2)))).'*exp(1i*(wave.omega*t))), ...
        'referenceCenter',referenceCenter,'referenceHeading',atan2(m.headingR(2,1),m.headingR(1,1)), ...
        'wave',wave,'waveFrequencyHz',wave.omega/(2*pi),'spectrum',spectrum,'parameters',p, ...
        'draft',h.frameDraftM,'equilibriumCenterOfMassZ',h.centerOfMassZ,'validity',validity, ...
        'hydroMode','6DOF composite glass/frame + equilibrium-waterline BEM + sidewall FK');
end

function p = validateParameters(p)
    positive = {'Lx','Ly','thickness','rho','g','duration','dt','fMinHz','fMaxHz','integrationStep'};
    for j = 1:numel(positive)
        validateattributes(p.(positive{j}),{'numeric'},{'scalar','real','finite','positive'}, ...
            'RIS_float',positive{j});
    end
    validateattributes(p.windSpeed,{'numeric'},{'scalar','real','finite','>=',3,'<=',25});
    validateattributes(p.windAngleDeg,{'numeric'},{'scalar','real','finite'});
    validateattributes(p.startTime,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(p.frequencyBins,{'numeric'},{'scalar','finite','integer','>=',2});
    validateattributes(p.directionBins,{'numeric'},{'scalar','finite','integer','>=',2});
    validateattributes(p.panelCount,{'numeric'},{'scalar','finite','integer','>=',4,'<=',80});
    validateattributes(p.radiationModes,{'numeric'},{'scalar','integer','>=',4});
    validateattributes(p.radiationFrequencyBins,{'numeric'},{'scalar','integer','>=',8});
    validatestring(p.initialization,{'waveSteady','rest'});
    sizes = [3,3,6]; names = {'initialPosition','initialEulerAngles','initialBodyVelocity'};
    for j = 1:3
        if ~isempty(p.(names{j}))
            validateattributes(p.(names{j}),{'numeric'},{'vector','numel',sizes(j),'real','finite'});
            p.(names{j}) = p.(names{j})(:);
        end
    end
    if mod(p.panelCount,2)~=0, error('RIS_float:Panels','panelCount must be even.'); end
    validateattributes(p.seed,{'numeric'},{'scalar','real','finite','integer','>=',0,'<=',2^32-1});
    pairs = {'center','horizontalInitialPosition','horizontalInitialVelocity'};
    for j = 1:numel(pairs)
        if j==3 && isempty(p.(pairs{j})), continue; end
        validateattributes(p.(pairs{j}),{'numeric'},{'vector','numel',2,'real','finite'});
        p.(pairs{j}) = p.(pairs{j})(:);
    end
    validateattributes(p.plot,{'logical','numeric'},{'scalar','real','finite'});
    if p.plot
        error('RIS_float:Plot','RIS_float only computes motion. Use run_RIS_wind_comparison for figures.');
    end
    if p.fMaxHz<=p.fMinHz, error('RIS_float:Band','fMaxHz must exceed fMinHz.'); end
end

function [w,s] = makeWave(p)
    if ~isempty(p.wave)
        w = p.wave;
        required = {'omega','kx','ky','coefficient'};
        if ~isstruct(w) || ~isscalar(w) || ~all(isfield(w,required))
            error('RIS_float:Wave','wave needs omega,kx,ky,coefficient vectors.');
        end
        for j = 1:4
            validateattributes(w.(required{j}),{'numeric'},{'vector','nonempty','finite'});
            w.(required{j}) = w.(required{j})(:);
        end
        validateattributes(w.omega,{'numeric'},{'real','positive'});
        validateattributes(w.kx,{'numeric'},{'real'});
        validateattributes(w.ky,{'numeric'},{'real'});
        if any([numel(w.kx),numel(w.ky),numel(w.coefficient)]~=numel(w.omega)) || ...
                any(abs(hypot(w.kx,w.ky)-w.omega.^2/p.g)>1e-9*max(1,w.omega.^2/p.g))
            error('RIS_float:Wave','Wave dimensions or deep-water dispersion are inconsistent.');
        end
        s = struct('name','supplied components','retainedFullSpectrumVariance',NaN);
    else
        [w,s] = RIS_wave_components(p);
    end
    if p.dt>pi/max(w.omega)
        error('RIS_float:TimeStep','dt exceeds the incident-wave Nyquist limit.');
    end
    s.waveVarianceM2 = sum(abs(w.coefficient).^2)/2;
    s.significantWaveHeightM = 4*sqrt(s.waveVarianceM2);
    s.rmsIncidentSlope = sqrt(sum(abs(w.coefficient).^2.*(w.omega.^2/p.g).^2)/2);
end

function [rao,A,B,F] = plateResponse(p,omega,kx,ky,index,M,K)
% Three parity sectors of Eq.(25); cache geometry/frequency, never sea phases.
    persistent keys entries
    if isempty(keys), keys = {}; entries = {}; end
    key = {p.Lx,p.Ly,p.innerLx,p.innerLy,p.rho,p.g,p.panelCount,omega};
    hit = find(cellfun(@(item) isequal(item,key),keys),1);
    if isempty(hit)
        [response,mass,damping,points] = solvePotentials(p,omega);
        if numel(keys)>=4, keys(1)=[]; entries(1)=[]; end
        keys{end+1} = key;
        entries{end+1} = struct('response',response,'mass',mass,'damping',damping,'points',points);
        hit = numel(keys);
    end
    entry = entries{hit};
    A = entry.mass;
    B = entry.damping;
    xy = entry.points;
    delta = p.Lx*p.Ly/p.panelCount^2;
    F = complex(zeros(3,numel(kx)));
    for j = 1:numel(omega)
        ids = find(index==j);
        if isempty(ids), continue; end
        cx = cos(xy(:,1)*kx(ids).'); sx = -1i*sin(xy(:,1)*kx(ids).');
        cy = cos(xy(:,2)*ky(ids).'); sy = -1i*sin(xy(:,2)*ky(ids).');
        r = entry.response(:,:,j);
        F(:,ids) = p.rho*p.g*4*delta*[r(:,1).'* (cx.*cy); ...
            r(:,2).'* (sx.*cy); r(:,3).'* (cx.*sy)];
    end
    Z = K-omega.'.^2.*(M+A)+1i*omega.'.*B;
    rao = F./Z(:,index);
end

function [R,A,B,xy] = solvePotentials(p,omega)
    n = p.panelCount;
    dx = p.Lx/n; dy = p.Ly/n;
    [ix,iy] = ndgrid(0:n/2-1,0:n/2-1);
    xy = [(ix(:)+0.5)*dx,(iy(:)+0.5)*dy];
    xy = xy(xy(:,1)>p.innerLx/2 | xy(:,2)>p.innerLy/2,:);
    [ox,oy] = ndgrid((0:n-1)*dx,(0:n-1)*dy);
    singular = zeros(n);
    for sx = [-1,1]
        for sy = [-1,1]
            xx = ox+sx*dx/2; yy = oy+sy*dy/2; rr = hypot(xx,yy);
            singular = singular+sx*sy*(xx.*log(yy+rr)+yy.*log(xx+rr))/(2*pi);
        end
    end
    [nodes,weights] = gaussLegendre(4);
    radius = zeros(n*n,16); quadrature = zeros(1,16);
    count = 0;
    for i = 1:4
        for j = 1:4
            count = count+1;
            radius(:,count) = reshape(hypot(ox+nodes(i)*dx/2,oy+nodes(j)*dy/2),[],1);
            quadrature(count) = dx*dy*weights(i)*weights(j)/4;
        end
    end
    mirrors = [1,1;-1,1;1,-1;-1,-1];
    ids = cell(1,4);
    for j = 1:4
        xx = round(abs(xy(:,1)-mirrors(j,1)*xy(:,1).')/dx)+1;
        yy = round(abs(xy(:,2)-mirrors(j,2)*xy(:,2).')/dy)+1;
        ids{j} = sub2ind([n,n],xx,yy);
    end
    P = [ones(size(xy,1),1),xy];
    parity = [1,1;-1,1;1,-1];
    R = complex(zeros(size(P,1),3,numel(omega)));
    A = zeros(3,numel(omega)); B = A;
    % H0(z)=2/pi*integral_0^(pi/2) sin(z*cos(theta)) dtheta.
    [gn,gw] = gaussLegendre(64);
    angles = (gn+1)*pi/4;
    for f = 1:numel(omega)
        k = omega(f)^2/p.g;
        z = k*radius;
        H0 = reshape(sin(z(:)*cos(angles).')*(gw/2),size(z));
        regular = -k/4*(bessely(0,z)+H0)-1i*k/2*besselj(0,z);
        table = singular+reshape(regular*quadrature.',n,n);
        for mode = 1:3
            sx = parity(mode,1); sy = parity(mode,2);
            H = table(ids{1})+sx*table(ids{2})+sy*table(ids{3})+sx*sy*table(ids{4});
            rr = (eye(size(H))+k*H)\P(:,mode);
            R(:,mode,f) = rr;
            T = p.rho*4*dx*dy*(P(:,mode).'*(P(:,mode)-rr))/k;
            A(mode,f) = real(T);
            B(mode,f) = -omega(f)*imag(T);
        end
    end
    if any(B(:)<-1e-6) || any(~isfinite(R(:)))
        error('RIS_float:Hydrodynamics','Nonpassive or nonfinite panel solution; refine panels/frequency band.');
    end
end

function [x,w] = gaussLegendre(n)
    b = (1:n-1)./sqrt(4*(1:n-1).^2-1);
    [V,D] = eig(diag(b,1)+diag(b,-1));
    [x,order] = sort(diag(D));
    w = 2*V(1,order).'.^2;
end
