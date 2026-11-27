function v = wris_validate(c,s)
%WRIS_VALIDATE Projection, complete-core orientation, clipping and convergence.
ns=c.validation.photons; reps=c.validation.repeats;
rows=zeros(numel(ns),4); powers=zeros(numel(ns),reps);
g=s.panels{1,2}.out.geometry;
for j=1:numel(ns)
    fprintf('Convergence: N=%d, repeats=%d\n',ns(j),reps);
    reported=zeros(1,reps);
    for k=1:reps
        cc=c; cc.mc.photons=ns(j); seed=c.validation.seed+100*j+k;
        out=optical(cc,g,seed);
        powers(j,k)=out.powerWithoutTurbulence; reported(k)=out.powerSE;
    end
    rows(j,:)=[ns(j),mean(powers(j,:)),std(powers(j,:)),mean(reported)];
end
v.convergence=array2table(rows,'VariableNames',{'photons','meanPowerW','powerStdW','reportedPowerSEW'});
v.powers=powers;
v.convergenceAgreement=abs(rows(end,2)-rows(end-1,2))< ...
    5*sqrt((rows(end,3)^2+rows(end-1,3)^2)/reps);
nAngles=size(s.panels,1);
v.transparent=cell(nAngles,1); ratioError=zeros(nAngles,1);
for row=1:nAngles
    cc=c; cc.mc.photons=max(ns); cc.water.absorption=0; cc.water.scattering=0;
    g=s.panels{row,1}.out.geometry;
    v.transparent{row}=optical(cc,g,901+row);
    ratioError(row)=abs(v.transparent{row}.directPlaneMetrics.ratio-1/cosd(g.thetaDeg));
end
frames=cellfun(@(p) p.out,s.panels,'UniformOutput',false);
v.maxEnergyResidual=max(cellfun(@(o) abs(o.energyResidual),frames),[],'all');
v.maxPowerIntegralErrorW=max(cellfun(@(o) ...
    abs(sum(o.noTurbulence,'all')*o.cellArea-o.powerWithoutTurbulence),frames),[],'all');
v.orientationErrorDeg=cellfun(@(o) abs(mod(o.directPlaneMetrics.orientationDeg- ...
    o.geometry.projectedOrientationDeg+90,180)-90),frames);
% Four sampling standard errors for an elliptical Gaussian principal axis.
v.orientationToleranceDeg=cellfun(@(o) 0.25+rad2deg(4/(sqrt(max(1, ...
    o.directPlanePowerWithoutTurbulence*o.photons/c.source.power))* ...
    max(o.directPlaneMetrics.ratio-1/o.directPlaneMetrics.ratio,0.02))),frames);
v.orientationDefined=cellfun(@(o) isfinite(o.geometry.projectedOrientationDeg) && ...
    isfinite(o.directPlaneMetrics.orientationDeg),frames);
v.maxProjectionRatioError=max(ratioError);
v.projectionRatioErrors=ratioError;
v.physicalAssertionsPassed=all(ratioError<0.05) && v.maxEnergyResidual<1e-10 && ...
    v.maxPowerIntegralErrorW<1e-10 && all(v.orientationErrorDeg(v.orientationDefined)<v.orientationToleranceDeg(v.orientationDefined)) && ...
    v.convergenceAgreement;
assert(v.physicalAssertionsPassed,'wris:Validation','Projection, rotation, energy or convergence check failed.');
writetable(v.convergence,fullfile(c.output,'convergence.csv'));
end

function out=optical(c,g,seed)
key=struct('source',c.source,'water',c.water,'mc',c.mc,'grid',c.grid, ...
    'size',[c.ris.Lx c.ris.Ly],'geometry',g,'seed',seed,'gain',1);
out=wris_cached(c,'optics',key,@() wris_mc(c,g,seed,1));
end
