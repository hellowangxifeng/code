function load = RIS_sidewall_load(p,wave,draft,zG)
%RIS_SIDEWALL_LOAD Exact first-order incident pressure integral on four walls.
% Output [Fx;Fy;Fz;Mx;My;Mz] per metre wave amplitude, body axes at equilibrium.
% p'=rho*g*exp(k*z)*exp(-i*(kx*x+ky*y)), -draft<=z<=0.
    arguments
        p (1,1) struct
        wave (1,1) struct
        draft (1,1) double {mustBePositive}
        zG (1,1) double {mustBeReal,mustBeFinite}
    end
    kx = wave.kx(:).'; ky = wave.ky(:).'; k = wave.omega(:).'.^2/p.g;
    [sx,dsx] = sincDerivative(kx*p.Lx/2);
    [sy,dsy] = sincDerivative(ky*p.Ly/2);
    area = p.Lx*p.Ly;
    footprint = area*sx.*sy;
    dx = area*p.Lx/2*dsx.*sy;
    dy = area*p.Ly/2*sx.*dsy;
    depth = -expm1(-k*draft)./k;
    % Weighted pressure depth relative to CG. Series avoids cancellation at kd=0.
    kd = k*draft;
    lever = -draft/2+zeros(size(k));
    large = kd>1e-3;
    lever(large) = 1./k(large)-draft./expm1(kd(large));
    lever(large) = -lever(large);
    lever(~large) = -draft/2+k(~large)*draft^2/12-k(~large).^3*draft^4/720;
    lever = lever-zG;
    factor = p.rho*p.g*depth;
    fx = 1i*factor.*kx.*footprint;
    fy = 1i*factor.*ky.*footprint;
    mz = factor.*(kx.*dy-ky.*dx);
    load = [fx;fy;zeros(size(fx));-lever.*fy;lever.*fx;mz];
end

function [s,ds] = sincDerivative(x)
    s = ones(size(x)); ds = zeros(size(x));
    small = abs(x)<1e-3;
    s(small) = 1-x(small).^2/6+x(small).^4/120;
    ds(small) = -x(small)/3+x(small).^3/30-x(small).^5/840;
    s(~small) = sin(x(~small))./x(~small);
    ds(~small) = (x(~small).*cos(x(~small))-sin(x(~small)))./x(~small).^2;
end
