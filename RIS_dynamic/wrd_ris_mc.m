function out = wrd_ris_mc(c,g,seed)
%WRD_RIS_MC Use W_RIS packets and existing quantized LC ray propagation.
[~,packets] = wris_mc(c,g,seed);
local = (packets.pos-g.pRis.')*g.RRis;
dirLocal = packets.dir*g.RRis;
qIn = c.optics.nWater*dirLocal(:,1:2);
kg = sqrt(c.optics.nGlass^2-sum(qIn.^2,2));
lc = local(:,1:2)+g.hIn*qIn./kg;
count = c.array.count;
pitch = c.array.pitchM;
ij = round(lc./pitch+(count+1)/2);
centres = (ij-(count+1)/2).*pitch;
use = all(ij>=1 & ij<=count,2) & all(abs(lc-centres)<=c.array.apertureM/2,2) & ...
    all(abs(lc)<=[c.ris.Lx c.ris.Ly]/2,2);
unitId = ij(use,1)+(ij(use,2)-1)*count(1);
[~,first,sampleUnit] = unique(unitId);
reference = centres(use,:);
reference = reference(first,:);
scale = c.source.power/c.mc.photons;
if isempty(reference)
    out = wrd_score(c,zeros(0,5));
    out.incidentPowerW = sum(packets.weight)*scale;
    out.transmittedPowerW = 0;
    out.activePowerW = 0;
    out.saturatedPowerFraction = NaN;
    out.energyResidual = packets.energyResidual;
    return
end

% Desired directions come from the physical receiver centre and the exact
% two-layer central rays. Apply the current electronics gradient bounds.
inPath = rp_layer_path(reference-g.sourceLocal(1:2).', ...
    [-g.sourceLocal(3),g.hIn],[c.optics.nWater,c.optics.nGlass]);
target = (c.receiver.center(:).'-g.pRis.')*g.RRis;
assert(target(3)>g.zTop,'wrd:Receiver','The RX centre must be above the RIS.');
outPath = rp_layer_path(target(1:2)-reference, ...
    [g.hOut,target(3)-g.zTop],[c.optics.nGlass,c.optics.nAir]);
[gradient,~,saturated] = rp_command(c,outPath.q-inPath.q);
phase = struct('referenceXY',reference,'gradientOverK0',gradient, ...
    'phaseCurvature',zeros(2,2,size(reference,1)),'phaseModel','linear');
weights = packets.weight(use).*rp_fresnel(c.optics.nWater,c.optics.nGlass,qIn(use,:));
illumination = struct('sampleUnitIndex',sampleUnit,'lcXY',lc(use,:), ...
    'qIn',qIn(use,:),'powerW',weights,'unitId',unitId);
trace = rp_trace(c,g,illumination,phase);
exit = trace.valid & trace.powerW>0;
ballistic = packets.ballistic(use);
hits = wrd_receive(c,trace.topWorld(exit,:),trace.directionWorld(exit,:), ...
    trace.powerW(exit),ballistic(exit));
out = wrd_score(c,hits);
out.energyResidual = packets.energyResidual;
out.incidentPowerW = sum(packets.weight)*scale;
out.activePowerW = sum(weights)*scale;
out.transmittedPowerW = sum(trace.powerW(exit))*scale;
out.saturatedPowerFraction = sum(weights(saturated(sampleUnit)))/max(sum(weights),realmin);
end
