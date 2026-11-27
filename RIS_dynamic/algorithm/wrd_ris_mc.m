function out = wrd_ris_mc(c,g,seed,sea)
%WRD_RIS_MC Moving-RIS link: W_RIS packets + existing quantized LC ray propagation.
%   out = wrd_ris_mc(c,g,seed)       no sea-surface diagnostics
%   out = wrd_ris_mc(c,g,seed,sea)   also reports RIS exit power that lies under
%                                    the undisturbed incident sea surface, and,
%                                    if c.comparison.risWaveBlocking is true,
%                                    removes RIS->RX rays that re-enter a crest.
%   Water transport is W_RIS/wris_mc unchanged. The LC tiles are re-pointed
%   every frame at the fixed RX (ideal per-frame update), with the device
%   gradient clamp and quantisation of RIS_phase/rp_command.
if nargin<4, sea = []; end
blocking = ~isempty(sea) && isfield(c.comparison,'risWaveBlocking') && c.comparison.risWaveBlocking;
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
incidentPowerW = sum(packets.weight)*scale;

if isempty(reference)
    out = finish(wrd_score(c,zeros(0,5)),packets.energyResidual,incidentPowerW,0,0,NaN,0,0);
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
overwashedW = 0;
blockedW = 0;
if ~isempty(sea) && any(exit)
    % The RIS optics assume air above the glass. Report how much exit power
    % leaves points that the undisturbed incident sea would cover.
    top = trace.topWorld(exit,:);
    power = trace.powerW(exit);
    overwashedW = sum(power(top(:,3)<sea.evaluate(top(:,1:2))))*scale;
end
if blocking
    [hits,blocked] = wrd_receive(c,trace.topWorld(exit,:),trace.directionWorld(exit,:), ...
        trace.powerW(exit),ballistic(exit),sea);
    blockedW = blocked*scale;
else
    hits = wrd_receive(c,trace.topWorld(exit,:),trace.directionWorld(exit,:), ...
        trace.powerW(exit),ballistic(exit));
end
out = finish(wrd_score(c,hits),packets.energyResidual,incidentPowerW, ...
    sum(weights)*scale,sum(trace.powerW(exit))*scale, ...
    sum(weights(saturated(sampleUnit)))/max(sum(weights),realmin),overwashedW,blockedW);
end

function out = finish(out,residual,incident,active,transmitted,saturated,overwashed,blocked)
% One fixed field order, so frames always concatenate into a struct array.
out.energyResidual = residual;
out.incidentPowerW = incident;
out.activePowerW = active;
out.transmittedPowerW = transmitted;
out.saturatedPowerFraction = saturated;
out.overwashedExitPowerW = overwashed;
out.waveBlockedPowerW = blocked;
end
