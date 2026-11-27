function h = RIS_flotation(p)
%RIS_FLOTATION Mass budget, exact upright displacement and composite inertia.
% Geometry/mass z coordinates are measured from the glass LOWER face.
% Payload and ballast are distributed symmetrically around the frame ring.
% frameDensity is an assumed effective dry density, including the frame skin.
    positive = {'Lx','Ly','thickness','glassDensity','frameWidth','frameHeight', ...
        'frameDensity','glassBottomAboveFrameBottom','rho','g'};
    for j = 1:numel(positive)
        validateattributes(p.(positive{j}),{'numeric'},{'scalar','real','finite','positive'});
    end
    validateattributes([p.payloadMass,p.ballastMass],{'numeric'},{'real','finite','nonnegative'});
    validateattributes([p.payloadZ,p.ballastZ],{'numeric'},{'real','finite'});
    mount = p.glassBottomAboveFrameBottom;
    if mount+p.thickness>p.frameHeight
        error('RIS_flotation:Mount','Glass must fit vertically inside the frame.');
    end
    outer = [p.Lx,p.Ly]+2*p.frameWidth;
    glassArea = p.Lx*p.Ly;
    frameArea = prod(outer)-glassArea;
    frameVolume = frameArea*p.frameHeight;
    masses = [glassArea*p.thickness*p.glassDensity; ...
        frameVolume*p.frameDensity;p.payloadMass;p.ballastMass];
    total = sum(masses);
    target = p.rho*frameArea*mount;
    % Continuous piecewise-linear Archimedes balance, including glass only
    % when the equilibrium waterline actually intersects it. No trapped air.
    if total<=target
        draft = total/(p.rho*frameArea);
    elseif total<=target+p.rho*(frameArea+glassArea)*p.thickness
        draft = mount+(total-target)/(p.rho*(frameArea+glassArea));
    else
        draft = (total/p.rho-glassArea*p.thickness)/frameArea;
    end
    if draft>=p.frameHeight
        error('RIS_flotation:Buoyancy','Total mass exceeds displacement before frame overtopping.');
    end
    glassDraft = max(0,min(p.thickness,draft-mount));
    lowerZ = mount-draft;
    z = [p.thickness/2;p.frameHeight/2-mount;p.payloadZ;p.ballastZ];
    cg = dot(masses,z)/total;
    % Ring second moments per unit area, exact outer minus inner rectangle.
    xx = (outer(1)*outer(2)^3-p.Lx*p.Ly^3)/(12*frameArea);
    yy = (outer(2)*outer(1)^3-p.Ly*p.Lx^3)/(12*frameArea);
    own = [p.Ly^2/12+p.thickness^2/12,p.Lx^2/12+p.thickness^2/12,(p.Lx^2+p.Ly^2)/12; ...
        xx+p.frameHeight^2/12,yy+p.frameHeight^2/12,xx+yy; ...
        xx,yy,xx+yy;xx,yy,xx+yy];
    own(:,1:2) = own(:,1:2)+(z-cg).^2;
    inertia = (masses.'*own).';
    volume = frameArea*draft+glassArea*glassDraft;
    zB = (frameArea*draft*(-draft/2)+ ...
        glassArea*glassDraft*(lowerZ+glassDraft/2))/volume;
    % At exact contact use the FRAME-only, dry-side tangent. The wet-side
    % stiffness is also reported; wet/dry switching is not solved dynamically.
    glassWet = glassDraft>1e-12 && glassDraft<p.thickness-1e-12;
    waterplaneArea = frameArea+glassWet*glassArea;
    waterplaneI = frameArea*[xx;yy]+glassWet*[p.Lx*p.Ly^3;p.Ly*p.Lx^3]/12;
    zG = lowerZ+cg;
    restoring = [p.rho*p.g*waterplaneArea; ...
        p.rho*p.g*waterplaneI+total*p.g*(zB-zG)];
    if any(restoring<=0), error('RIS_flotation:Stability','Nonpositive hydrostatic stiffness.'); end
    h = struct('outerSizeM',outer,'glassAreaM2',glassArea,'frameWaterplaneAreaM2',frameArea, ...
        'frameVolumeM3',frameVolume,'glassMassKg',masses(1),'frameMassKg',masses(2), ...
        'payloadMassKg',masses(3),'ballastMassKg',masses(4),'totalMassKg',total, ...
        'targetContactMassKg',target,'requiredBallastForContactKg',target-sum(masses(1:3)), ...
        'topAtWaterlineMassKg',target+p.rho*prod(outer)*p.thickness, ...
        'frameOnlyMaxDisplacementKg',p.rho*frameVolume, ...
        'frameDraftM',draft,'frameFreeboardM',p.frameHeight-draft, ...
        'glassDraftM',glassDraft,'glassBottomZ',lowerZ,'glassTopZ',lowerZ+p.thickness, ...
        'glassBottomFromCGM',-cg,'glassTopFromCGM',p.thickness-cg, ...
        'centerOfMassZ',zG,'centerOfBuoyancyZ',zB,'inertiaKgM2',inertia, ...
        'restoring',restoring,'frameWaterplaneSecondMomentM4',frameArea*[xx;yy], ...
        'contactHeaveStiffnessDryWetNpm',p.rho*p.g*[frameArea,prod(outer)], ...
        'buoyancyResidualN',p.rho*p.g*volume-total*p.g, ...
        'bottomTouchesCalmWater',abs(lowerZ)<1e-10,'topInAir',lowerZ+p.thickness>0);
end
