function I=rp_direct_density(c,g,worldXY)
%RP_DIRECT_DENSITY Exact ray-area Jacobian for continuous Gaussian direct beam.
depth=-g.sourceLocal(3); height=g.rxPlaneZ-g.pRis(3);
target=worldXY-g.pRis(1:2).';
path=rp_layer_path(target-g.sourceLocal(1:2).',[depth height], ...
    [c.optics.nWater,c.optics.nAir],c.numerics.bisectionIterations);
q=path.q; kw=path.kz(:,1); ka=path.kz(:,2); q2=sum(q.^2,2);
bottom=target-height*q./ka;
r=([bottom,zeros(size(bottom,1),1)]-g.sourceLocal.')*g.RTx;
axial=r(:,3); lengthWater=vecnorm(r,2,2); normal=g.RTx.'*[0;0;1];
I=2*c.source.power/(pi*c.source.divergence^2)* ...
    exp(-2*sum(r(:,1:2).^2,2)./(c.source.divergence^2*axial.^2) ...
    -(c.water.absorption+c.water.scattering)*lengthWater).*max(r*normal,0)./axial.^3;
aa=depth./kw+height./ka; bb=depth./kw.^3+height./ka.^3;
detWater=(depth./kw).^2.*(1+q2./kw.^2);
I=I.*detWater./(aa.*(aa+bb.*q2)).*rp_fresnel(c.optics.nWater,c.optics.nAir,q);
I(axial<=0)=0;
end
