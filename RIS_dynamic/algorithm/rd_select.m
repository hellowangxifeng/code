function choice = rd_select(A,O,M,family,feasible,policy,design)
%RD_SELECT Design chosen by each scheme at each transmit power (RIS_phase rules).
%   A, O, M: designs x powers ABER, OP and mean received power.
%   choice(scheme, power, metric), metric 1 = ABER figure, 2 = OP figure;
%   scheme 1 = Direct (0, no design), 2 Steering-only, 3 Flat top,
%   4 Robust ABER (metric 1) / Robust OP (metric 2), 5 Max mean power.
%     Flat top       min ABER (OP) over feasible flat_top designs
%     Robust         min ABER (OP) over feasible robust designs
%     Max mean power max <P_r> over feasible flat_top and robust designs
%   Ties (e.g. OP = 1 everywhere) are broken by the other risk, then by the
%   larger mean power (RIS_phase rp_fixed_receiver_collect/choosePhase).
%   policy 'design_power' freezes every choice at the power index design.
np = size(A,2);
steer = find(family=="steering",1);
flat = pool(family=="flat_top",feasible);
robust = pool(family=="robust",feasible);
shaped = pool(ismember(family,["flat_top","robust"]),feasible);
choice = zeros(5,np,2);
choice(2,:,:) = steer;
choice(3,:,1) = best(A,O,M,flat);
choice(3,:,2) = best(O,A,M,flat);
choice(4,:,1) = best(A,O,M,robust);
choice(4,:,2) = best(O,A,M,robust);
score = M; score(~shaped,:) = -Inf;
[~,k] = max(score,[],1);
choice(5,:,1) = k; choice(5,:,2) = k;
if strcmp(policy,'design_power')
    choice = repmat(choice(:,design,:),1,np,1);
end
end

function allowed = pool(member,feasible)
allowed = member(:) & feasible(:);
if ~any(allowed), allowed = member(:); end      % report infeasibility, never crash
assert(any(allowed),'rd:Pool','A design family has no candidate.');
end

function ids = best(primary,secondary,power,allowed)
ids = zeros(1,size(primary,2));
for j = 1:numel(ids)
    p = find(allowed);
    b = min(primary(p,j)); p = p(primary(p,j)<=b+8*eps(max(1,abs(b))));
    b = min(secondary(p,j)); p = p(secondary(p,j)<=b+8*eps(max(1,abs(b))));
    [~,k] = max(power(p,j)); ids(j) = p(k);
end
end
