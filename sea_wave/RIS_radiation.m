function radiation = RIS_radiation(omega,addedMass,damping,modeCount)
%RIS_RADIATION Passive rational approximation of BEM radiation impedance.
% Fits B(w)+i*w*A(w) with i*w*Ainf + sum c*(s+a)/((s+a)^2+b^2).
% Nonnegative residues and Ainf guarantee passivity; errors are reported.
    arguments
        omega (:,1) double {mustBePositive,mustBeFinite}
        addedMass (3,:) double {mustBeFinite}
        damping (3,:) double {mustBeFinite}
        modeCount (1,1) double {mustBeInteger,mustBeGreaterThanOrEqual(modeCount,4)}
    end
    frequency = logspace(log10(min(omega)/2),log10(max(omega)*2),modeCount).';
    % Pole spacing/width is a numerical approximation choice, not a drag coefficient.
    width = frequency*(exp(log(frequency(end)/frequency(1))/(modeCount-1))-1);
    s = 1i*omega;
    basis = (s+width.')./((s+width.').^2+frequency.'.^2);
    design = [basis,s];
    residues = zeros(modeCount,3); infinity = zeros(3,1);
    fittedMass = zeros(size(addedMass)); fittedDamping = fittedMass;
    errorRelative = zeros(3,1);
    for j = 1:3
        target = damping(j,:).'+s.*addedMass(j,:).';
        scale = max(abs(target),0.02*max(abs(target)));
        weighted = design./scale;
        realDesign = [real(weighted);imag(weighted)];
        columnScale = vecnorm(realDesign);
        solution = lsqnonneg(realDesign./columnScale,[real(target./scale);imag(target./scale)]);
        solution = solution./columnScale.';
        residues(:,j) = solution(1:end-1);
        infinity(j) = solution(end);
        fit = design*solution;
        fittedMass(j,:) = imag(fit)./omega;
        fittedDamping(j,:) = real(fit);
        errorRelative(j) = norm((fit-target)./scale)/sqrt(numel(omega));
    end
    radiation = struct('frequency',frequency,'width',width,'residues',residues, ...
        'infiniteAddedMass',infinity,'fittedAddedMass',fittedMass, ...
        'fittedDamping',fittedDamping,'relativeImpedanceRmsError',errorRelative, ...
        'passive',all(residues(:)>=0) && all(infinity>=0));
end
