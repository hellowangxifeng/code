function distance = rd_crossing(p,d,limit,sea,settings,sense)
%RD_CROSSING First water/air crossing on each finite straight segment.
% sense=+1 searches from water into air, -1 from air into water.
% Bracket the possible boundary slab. Search nonmonotone rays in path order.
distance = inf(size(limit));
if sea.amplitudeBound == 0
    s = -p(:,3)./d(:,3);
    hit = sense*d(:,3)>0 & s>=0 & s<=limit;
    distance(hit) = s(hit);
    return
end

dz = d(:,3);
horizontal = dz==0;
edge = sea.amplitudeBound;
s1 = (-edge-p(:,3))./dz;
s2 = ( edge-p(:,3))./dz;
start = max(0,min(s1,s2));
stop = min(limit,max(s1,s2));
start(horizontal) = 0;
stop(horizontal) = limit(horizontal);
slope = abs(d(:,1:2))*sea.slopeBound.';
possible = stop>=start & stop>0 & sense*dz+slope>0 & ...
    (~horizontal | abs(p(:,3))<=edge);
if ~isempty(sea.gridBounds)
    possibleIds = find(possible);
    a = p(possibleIds,1:2)+reshape(start(possibleIds),[],1).*d(possibleIds,1:2);
    b = p(possibleIds,1:2)+reshape(stop(possibleIds),[],1).*d(possibleIds,1:2);
    lower = sea.gridBounds(:,1).';
    upper = sea.gridBounds(:,2).';
    onGrid = all(min(a,b)>=lower & max(a,b)<=upper,2);
    gridIds = possibleIds(onGrid);
    slope(gridIds) = abs(d(gridIds,1:2))*sea.gridSlopeBound.';
    nonHorizontal = gridIds(~horizontal(gridIds));
    localEdge = sea.gridAmplitudeBound;
    s1 = (-localEdge-p(nonHorizontal,3))./dz(nonHorizontal);
    s2 = ( localEdge-p(nonHorizontal,3))./dz(nonHorizontal);
    start(nonHorizontal) = max(start(nonHorizontal),min(s1,s2));
    stop(nonHorizontal) = min(stop(nonHorizontal),max(s1,s2));
    possible(gridIds) = stop(gridIds)>=start(gridIds) & sense*dz(gridIds)+slope(gridIds)>0 ...
        & (~horizontal(gridIds) | abs(p(gridIds,3))<=localEdge);
end
ids = find(possible);
if isempty(ids), return; end

left = start(ids);
right = stop(ids);
gap = signedGap(p(ids,:),d(ids,:),left,sea,sense);
first = gap>=0;
left(first) = max(0,left(first)-settings.rootTolerance);
right(first) = start(ids(first));
monotone = sense*dz(ids)>slope(ids);
checkEnd = find(monotone & ~first);
first(checkEnd) = signedGap(p(ids(checkEnd),:),d(ids(checkEnd),:), ...
    right(checkEnd),sea,sense)>=0;

remaining = find(~monotone & ~first);
while ~isempty(remaining)
    % A derivative bound allows larger steps far from the interface. The
    % minimum step resolves grazing crossings within the stated resolution.
    ray = ids(remaining);
    bound = sense*dz(ray)+slope(ray);
    step = max(settings.marchStep, -gap(remaining)./bound);
    next = min(stop(ray),left(remaining)+step);
    value = signedGap(p(ray,:),d(ray,:),next,sea,sense);
    crossed = value>=0;
    right(remaining(crossed)) = next(crossed);
    first(remaining(crossed)) = true;
    unfinished = ~crossed & next<stop(ray);
    remaining = remaining(unfinished);
    left(remaining) = next(unfinished);
    gap(remaining) = value(unfinished);
end

ids = ids(first);
left = left(first);
right = right(first);
while any(right-left>settings.rootTolerance)
    mid = (left+right)/2;
    onTargetSide = signedGap(p(ids,:),d(ids,:),mid,sea,sense)>=0;
    right(onTargetSide) = mid(onTargetSide);
    left(~onTargetSide) = mid(~onTargetSide);
end
distance(ids) = right;
end

function value = signedGap(p,d,s,sea,sense)
xyz = p+d.*s(:);
value = sense*(xyz(:,3)-sea.evaluate(xyz(:,1:2)));
end
