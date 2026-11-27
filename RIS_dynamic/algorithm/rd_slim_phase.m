function s = rd_slim_phase(p,family,width,margin)
%RD_SLIM_PHASE Keep only what rp_trace needs from a RIS_phase phase design.
% The LC commands are frozen in the BOARD frame; under waves the same
% per-tile gradients are applied to the moving board.
n = numel(p.unitIds);
s = struct('unitIds',p.unitIds(:),'referenceXY',p.referenceXY, ...
    'gradientOverK0',p.gradientOverK0,'phaseCurvature',zeros(2,2,n), ...
    'phaseModel','linear','commandDeg',p.commandDeg,'family',family, ...
    'targetWidthM',width,'margin',margin);
end
