function [xy,valid] = rp_carrier_landing(c,g,reference,q)
%RP_CARRIER_LANDING Exact layered ray from one phase-screen point to world RX plane.
arguments
    c (1,1) struct
    g (1,1) struct
    reference (1,2) double {mustBeReal,mustBeFinite}
    q (:,2) double {mustBeReal,mustBeFinite}
end
air=c.optics.nAir^2-sum(q.^2,2);
valid=air>0;
xy=nan(size(q));
if ~any(valid), return; end
ids=find(valid); selected=q(ids,:);
glass=sqrt(c.optics.nGlass^2-sum(selected.^2,2));
top=[reference+g.hOut*selected./glass,repmat(g.zTop,numel(ids),1)];
world=top*g.RRis.'+g.pRis.';
direction=[selected,sqrt(air(ids))]/c.optics.nAir*g.RRis.';
distance=(g.rxPlaneZ-world(:,3))./max(direction(:,3),realmin);
forward=direction(:,3)>c.numerics.forwardCosTolerance & distance>0;
landing=world+distance.*direction;
xy(ids(forward),:)=landing(forward,1:2);
valid(ids(~forward))=false;
end
