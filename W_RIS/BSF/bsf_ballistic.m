function intensity = bsf_ballistic(c, g, x, y)
%BSF_BALLISTIC Exact zero-collision ray density for the ORIGINAL Gaussian source.
% Finite-waist expectation is integrated by a 3x3 Gaussian cubature.
% x,y are local RIS coordinates with matching sizes, not coordinate vectors.
arguments
    c (1,1) struct
    g (1,1) struct
    x {mustBeNumeric, mustBeReal, mustBeFinite}
    y {mustBeNumeric, mustBeReal, mustBeFinite}
end
assert(isequal(size(x),size(y)), 'bsf:Shape', 'x and y must have matching sizes.');
q = [x(:),y(:),zeros(numel(x),1)]*g.RRis.' + g.pRis.';
r = (q-g.pTx.')*g.RTx;
normal = g.RTx.'*g.normal;
z = r(:,3);
valid = z>0;
intensity = zeros(numel(x),1);
% N(0,1): nodes {0,+/-sqrt(3)}, weights {2/3,1/6,1/6}.
nodes = [0,sqrt(3),-sqrt(3)]*c.source.waist/2;
weights = [2/3,1/6,1/6];
sigmaAngle = c.source.divergence/2;
extinction = c.water.absorption+c.water.scattering;
for i = 1:3
    for j = 1:3
        rr = r(valid,:) - [nodes(i),nodes(j),0];
        zz = z(valid);
        rayLength = vecnorm(rr,2,2);
        % d(alpha_x,alpha_y)/dA = (n dot rr)/z^3, NOT cos(mu)/z^2.
        jacobian = max(rr*normal,0)./zz.^3;
        density = exp(-sum(rr(:,1:2).^2,2)./(2*sigmaAngle^2*zz.^2))/(2*pi*sigmaAngle^2);
        density(rayLength>c.mc.maxPath) = 0;
        intensity(valid) = intensity(valid) + weights(i)*weights(j)*c.source.power* ...
            density.*jacobian.*exp(-extinction*rayLength);
    end
end
intensity = reshape(intensity,size(x));
end
