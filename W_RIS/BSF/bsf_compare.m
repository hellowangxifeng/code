function metrics = bsf_compare(prediction, reference, x, y, area, modelSE, referenceSE)
%BSF_COMPARE Absolute-power, spatial, and noise-aware comparisons on SAME bins.
arguments
    prediction {mustBeNumeric, mustBeReal, mustBeFinite}
    reference {mustBeNumeric, mustBeReal, mustBeFinite}
    x {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
    y {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
    area (1,1) {mustBeNumeric, mustBeReal, mustBePositive}
    modelSE {mustBeNumeric, mustBeReal, mustBeFinite}
    referenceSE {mustBeNumeric, mustBeReal, mustBeFinite}
end
assert(isequal(size(prediction),size(reference),size(modelSE),size(referenceSE)), ...
    'bsf:Shape', 'All intensity and SE maps must have matching sizes.');
[X,Y] = ndgrid(x,y);
difference = prediction-reference;
p = sum(prediction,'all')*area;
r = sum(reference,'all')*area;
metrics.powerW = p;
metrics.referencePowerW = r;
metrics.powerRelative = (p-r)/r;
metrics.relativeL2 = norm(difference,'fro')/norm(reference,'fro');
metrics.relativeL1 = sum(abs(difference),'all')/sum(reference,'all');
metrics.peakRelative = (max(prediction,[],'all')-max(reference,[],'all'))/max(reference,[],'all');
metrics.maxAbsoluteOverPeak = max(abs(difference),[],'all')/max(reference,[],'all');
metrics.referenceNoiseL2 = norm(referenceSE,'fro')/norm(reference,'fro');
metrics.modelNoiseL2 = norm(modelSE,'fro')/norm(prediction,'fro');
metrics.amplitudeL2 = norm(sqrt(max(prediction,0))-sqrt(max(reference,0)),'fro')/ ...
    norm(sqrt(max(reference,0)),'fro');
muP = [sum(X.*prediction,'all'),sum(Y.*prediction,'all')]/sum(prediction,'all');
muR = [sum(X.*reference,'all'),sum(Y.*reference,'all')]/sum(reference,'all');
metrics.centroidM = muP;
metrics.referenceCentroidM = muR;
metrics.centroidErrorM = norm(muP-muR);
metrics.covarianceM2 = covariance2(X,Y,prediction,muP);
metrics.referenceCovarianceM2 = covariance2(X,Y,reference,muR);
% Normal-score diagnostics only where the reference effective count >=25.
variance = modelSE.^2+referenceSE.^2;
mask = reference>=5*referenceSE & referenceSE>0 & variance>0;
metrics.eligibleBins = nnz(mask);
metrics.meanSquaredZ = mean(difference(mask).^2./variance(mask));
metrics.coverage95 = mean(abs(difference(mask))<=1.96*sqrt(variance(mask)));
metrics.minimumIntensity = min(prediction,[],'all');
end

function covariance = covariance2(X,Y,I,mu)
mass = sum(I,'all');
xx = sum((X-mu(1)).^2.*I,'all')/mass;
yy = sum((Y-mu(2)).^2.*I,'all')/mass;
xy = sum((X-mu(1)).*(Y-mu(2)).*I,'all')/mass;
covariance = [xx,xy;xy,yy];
end
