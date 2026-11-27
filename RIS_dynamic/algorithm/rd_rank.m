function r = rd_rank(curves,powerDBm,c)
%RD_RANK Required transmit power per scheme and the clearly best RIS optimizer(s).
%   For metric ABER (OP) the required power P_req(s) is the first P_t at
%   which scheme s reaches c.rank.targets(1) (2), by log-linear
%   interpolation; Inf if it never does inside c.metrics.powerDBm.
%   Among the three optimizers (3 Flat top, 4 Robust, 5 Max mean power) the
%   winners are those within c.rank.toleranceDB of the best one. clear is
%   true when every winner is at least c.rank.marginDB ahead of EVERY other
%   link (the remaining optimizers, Steering-only and Direct). If no
%   optimizer reaches the target, the metric at the highest power decides.
names = {'aber','op'};
for k = 1:2
    y = curves.(names{k});
    need = inf(1,5);
    for s = 1:5
        need(s) = crossing(powerDBm,y(s,:),c.rank.targets(k));
    end
    optimizers = 3:5;
    [bestNeed,~] = min(need(optimizers));
    if isfinite(bestNeed)
        winners = optimizers(need(optimizers)<=bestNeed+c.rank.toleranceDB);
        others = setdiff(1:5,winners);
        margin = min(need(others))-max(need(winners));
    else
        last = log10(max(y(optimizers,end),realmin));
        winners = optimizers(last<=min(last)+1e-9);
        margin = 0;
    end
    r.(names{k}) = struct('requiredDBm',need,'winners',winners,'marginDB',margin, ...
        'clear',isfinite(bestNeed) && margin>=c.rank.marginDB, ...
        'gainOverSteeringDB',need(2)-bestNeed,'gainOverDirectDB',need(1)-bestNeed, ...
        'target',c.rank.targets(k));
end
end

function P = crossing(powerDBm,y,target)
P = Inf;
y = double(y(:)).';
k = find(y<=target,1);
if isempty(k), return; end
if k==1, P = powerDBm(1); return; end
a = log10(max(y(k-1),realmin)); b = log10(max(y(k),realmin));
t = (log10(target)-a)/(b-a);
P = powerDBm(k-1)+t*(powerDBm(k)-powerDBm(k-1));
end
