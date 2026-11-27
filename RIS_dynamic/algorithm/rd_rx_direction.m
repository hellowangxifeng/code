function q0 = rd_rx_direction(c,a,base)
%RD_RX_DIRECTION Board-local exit direction of the common receiver (direction cosines).
%   Every tile m can only emit q = q_in,m + b with |b_x|,|b_y| <= B
%   (B = sin 15.414 deg, rp_command). Per axis, the directions reachable by
%   all tiles except a 1e-5 power tail form [l_d, h_d] (RIS_phase eq. 47):
%     l_d = Q_{1-eps}(q_in,d - B),  h_d = Q_{eps}(q_in,d + B).
%   q0 is the point of that box closest to the board normal, moved
%   c.rx.carrierMargin further inside so that finite-width shaping and
%   jitter-robust spreading keep their tiles unsaturated. The RX is placed
%   where the ray from the LC centre with this q0 meets the RX plane.
ids = base.unitIds;
weights = a.unitDirectPowerW(ids);
weights = weights/max(sum(weights),realmin);
bound = c.optics.nAir*sind(floor(c.device.maxAngleDeg/c.device.angleStepDeg)*c.device.angleStepDeg);
tail = c.study.maximumTailPower/4;
q0 = zeros(1,2);
for d = 1:2
    lo = weightedQuantile(base.qInAtReference(:,d)-bound,weights,1-tail);
    hi = weightedQuantile(base.qInAtReference(:,d)+bound,weights,tail);
    assert(lo<=hi,'rd:Reach','No direction is reachable by all illuminated tiles.');
    q0(d) = min(hi,max(lo,0));
    step = min(c.rx.carrierMargin,(hi-lo)/2);
    if q0(d)==lo && lo>0
        q0(d) = q0(d)+step;
    elseif q0(d)==hi && hi<0
        q0(d) = q0(d)-step;
    end
end
assert(sum(q0.^2)<c.optics.nAir^2,'rd:Reach','The receiver direction is evanescent in air.');
end

function v = weightedQuantile(x,w,f)
[x,order] = sort(x);
t = cumsum(w(order));
v = x(find(t>=f*t(end),1));
end
