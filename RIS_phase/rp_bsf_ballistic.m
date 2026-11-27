function intensity = rp_bsf_ballistic(c,g,xy,height)
%RP_BSF_BALLISTIC PDF Eq.(8), including finite waist and exact ray-area mapping.
% Bottom density is integrated over the 3x3 Gaussian source cubature. Inside
% glass each source node follows its own Snell path and area Jacobian.
arguments
    c (1,1) struct
    g (1,1) struct
    xy (:,2) double {mustBeReal,mustBeFinite}
    height (1,1) double {mustBeNonnegative} = 0
end
nw=c.optics.nWater; ng=c.optics.nGlass;
sourceRotation=g.RRis.'*g.RTx;
nodes=[0 sqrt(3) -sqrt(3)]*c.source.waist/2;
weights=[2/3 1/6 1/6];
intensity=zeros(size(xy,1),1);
for i=1:3
    for j=1:3
        origin=g.sourceLocal+sourceRotation*[nodes(i);nodes(j);0];
        depth=-origin(3);
        path=rp_layer_path(xy-origin(1:2).',[depth height],[nw ng],c.numerics.bisectionIterations);
        q=path.q; kw=path.kz(:,1); kg=path.kz(:,2);
        bottom=xy-height*q./kg;
        r=([bottom,zeros(size(xy,1),1)]-origin.')*sourceRotation;
        z=r(:,3); ell=vecnorm(r,2,2);
        normal=sourceRotation.'*[0;0;1];
        good=z>0 & ell<=c.mc.maxPath & all(abs(bottom)<=[c.ris.Lx c.ris.Ly]/2,2);
        if height==0
            % The lower-plane density is defined outside the finite board too.
            good=z>0 & ell<=c.mc.maxPath;
        end
        density=zeros(size(z));
        density(good)=2/(pi*c.source.divergence^2)* ...
            exp(-2*sum(r(good,1:2).^2,2)./(c.source.divergence^2*z(good).^2) ...
            -(c.water.absorption+c.water.scattering)*ell(good)).* ...
            max(r(good,:)*normal,0)./z(good).^3;
        if height>0
            q2=sum(q.^2,2);
            detWater=(depth./kw).^2.*(1+q2./kw.^2);
            aa=depth./kw+height./kg;
            bb=depth./kw.^3+height./kg.^3;
            density=density.*detWater./(aa.*(aa+bb.*q2)).*rp_fresnel(nw,ng,q);
        end
        intensity=intensity+c.source.power*weights(i)*weights(j)*density;
    end
end
end
