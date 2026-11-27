function [gradient,commandDeg,saturated] = rp_command(c,ideal)
%RP_COMMAND Apply the same independent-axis clamp and quantisation to all designs.
arguments
    c (1,1) struct
    ideal (:,2) double {mustBeReal,mustBeFinite}
end
bound=c.optics.nAir*sind(c.device.maxAngleDeg);
projected=max(-bound,min(bound,ideal));
limit=floor(c.device.maxAngleDeg/c.device.angleStepDeg+1e-12);
commandDeg=max(-limit,min(limit,round(asind(projected/c.optics.nAir)/c.device.angleStepDeg))) ...
    *c.device.angleStepDeg;
gradient=c.optics.nAir*sind(commandDeg);
% Beyond the largest implemented grid point, even if below the continuous limit.
effectiveBound=c.optics.nAir*sind(limit*c.device.angleStepDeg);
saturated=any(abs(ideal)>effectiveBound+1e-12,2);
end
