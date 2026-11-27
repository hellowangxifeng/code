function [hits,blockedWeight] = wrd_receive(c,p,d,weights,ballistic,sea)
%WRD_RECEIVE Intersect the same physical, tilted detector plane in both links.
normal = c.receiver.normal(:)/norm(c.receiver.normal);
if abs(normal(3))<0.9
    u = cross([0;0;1],normal);
else
    u = cross([0;1;0],normal);
end
u = u/norm(u);
v = cross(normal,u);
cosine = d*normal;
distance = ((c.receiver.center(:).'-p)*normal)./cosine;
forward = cosine>1e-10 & distance>0;
ids = find(forward);
delta = p(ids,:)+distance(ids).*d(ids,:)-c.receiver.center(:).';
xy = delta*[u,v];
edges = c.receiver.mapEdges;
scored = all(xy>=edges(1) & xy<=edges(end),2) | sum(xy.^2,2)<=c.receiver.radius^2;
ids = ids(scored);
xy = xy(scored,:);
blockedWeight = 0;
if nargin>=6 && ~isempty(ids)
    % Start just inside air to avoid counting the transmitted interface twice.
    offset = 8*c.surface.rootTolerance;
    search = max(0,distance(ids)-offset);
    recross = wrd_crossing(p(ids,:)+offset*d(ids,:),d(ids,:),search,sea,c.surface,-1);
    blocked = isfinite(recross);
    blockedWeight = sum(weights(ids(blocked)));
    ids = ids(~blocked);
    xy = xy(~blocked,:);
end
hits = [xy,weights(ids),double(ballistic(ids)),cosine(ids)];
end
