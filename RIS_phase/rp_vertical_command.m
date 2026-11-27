function [gradient,commandDeg,directionWorld,valid] = rp_vertical_command(c,g,qIn)
%RP_VERTICAL_COMMAND Globally closest world +z carrier on the discrete command grid.
% Enumerate X. For each X the concave objective has one Y maximum; its two
% adjacent realizable Y commands suffice. Complexity is linear in axis count.
arguments
    c (1,1) struct
    g (1,1) struct
    qIn (:,2) double {mustBeReal,mustBeFinite}
end
vertical = g.RRis(3,:);
assert(vertical(3)>0,'rp:Pose','The RIS upper normal must have a positive world z component.');
nAir = c.optics.nAir;
limit = floor(c.device.maxAngleDeg/c.device.angleStepDeg+1e-12);
commandIndex = (-limit:limit).';
bAxis = nAir*sind(commandIndex*c.device.angleStepDeg);
[gradient,commandDeg] = rp_command(c,-qIn); % legal fallback if every mode is evanescent
directionWorld = nan(size(qIn,1),3);
valid = false(size(qIn,1),1);
for m = 1:size(qIn,1)
    qx = qIn(m,1)+bAxis;
    disk = nAir^2-qx.^2;
    use = disk>0;
    if ~any(use), continue; end
    qx = qx(use); xIndex = commandIndex(use); radius = sqrt(disk(use));
    % Restrict the Y grid to propagating modes, excluding the grazing boundary.
    lower = max(-1,min(1,(-radius-qIn(m,2))/nAir));
    upper = max(-1,min(1,( radius-qIn(m,2))/nAir));
    lo = max(-limit,ceil(asind(lower)/c.device.angleStepDeg+1e-10));
    hi = min( limit,floor(asind(upper)/c.device.angleStepDeg-1e-10));
    bestY = vertical(2)*radius/hypot(vertical(2),vertical(3));
    continuous = asind(max(-1,min(1,(bestY-qIn(m,2))/nAir)))/c.device.angleStepDeg;
    yIndex = max(lo,min(hi,[floor(continuous),ceil(continuous)]));
    qy = qIn(m,2)+nAir*sind(yIndex*c.device.angleStepDeg);
    kz = sqrt(max(0,nAir^2-qx.^2-qy.^2));
    score = (vertical(1)*qx+vertical(2)*qy+vertical(3)*kz)/nAir;
    legal = lo<=hi & nAir^2-qx.^2-qy.^2>0;
    score(~legal | score<=c.numerics.forwardCosTolerance) = -Inf;
    [value,index] = max(score,[],'all','linear');
    if ~isfinite(value), continue; end
    [row,column] = ind2sub(size(score),index);
    commandDeg(m,:) = [xIndex(row),yIndex(row,column)]*c.device.angleStepDeg;
    gradient(m,:) = nAir*sind(commandDeg(m,:));
    q = qIn(m,:)+gradient(m,:);
    directionWorld(m,:) = [q,sqrt(nAir^2-sum(q.^2))]/nAir*g.RRis.';
    valid(m) = true;
end
end
