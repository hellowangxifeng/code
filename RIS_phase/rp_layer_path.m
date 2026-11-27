function r = rp_layer_path(deltaXY,heights,indices,iterations)
%RP_LAYER_PATH Exact parallel-layer Fermat ray with variable endpoint heights.
arguments
    deltaXY (:,2) double {mustBeReal,mustBeFinite}
    heights (:,:) double {mustBeReal,mustBeNonnegative,mustBeFinite}
    indices (1,:) double {mustBeReal,mustBePositive,mustBeFinite}
    iterations (1,1) double {mustBeReal,mustBeInteger,mustBePositive} = 48
end
assert(size(heights,2)==numel(indices) && all(any(heights>0,2)),'rp:Layers','Invalid layers.');
assert(size(heights,1)==1 || size(heights,1)==size(deltaXY,1),'rp:Layers','Wrong layer count.');
if size(heights,1)==1, heights = repmat(heights,size(deltaXY,1),1); end
rho = vecnorm(deltaXY,2,2);
unit = deltaXY./max(rho,realmin);
lo = zeros(size(rho));
positiveIndex = repmat(indices,size(deltaXY,1),1);
positiveIndex(heights<=0) = Inf;
hi = min(positiveIndex,[],2)*(1-1e-13);
for j = 1:iterations
    q = (lo+hi)/2;
    kz = sqrt(max(indices.^2-q.^2,realmin));
    distance = sum(heights.*q./kz,2);
    lower = distance<rho;
    lo(lower) = q(lower); hi(~lower) = q(~lower);
end
q = (lo+hi)/2; q(rho==0) = 0;
kz = sqrt(max(indices.^2-q.^2,realmin));
r.q = unit.*q;
r.kz = kz;
r.length = heights.*indices./kz;
r.opl = sum(r.length.*indices,2);
r.residual = sum(heights.*q./kz,2)-rho;
end
