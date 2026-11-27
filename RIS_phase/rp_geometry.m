function g = rp_geometry(c,m,tx,rx,k)
%RP_GEOMETRY Freeze pose and finite-source geometry; no chief-ray aperture gate.
arguments
    c (1,1) struct
    m (1,1) struct
    tx (1,1) struct
    rx (1,1) struct
    k (1,1) {mustBeNumeric,mustBeReal,mustBeInteger,mustBePositive}
end
g.RRis = m.rotationMatrix(:,:,k); g.pRis = m.bottomCenter(:,k);
g.pTx = tx.position(:,k); g.dTx = tx.direction(:,k);
g.RTx = tx.rotation(:,:,k); g.pBottom = g.pRis; g.normal = g.RRis(:,3);
g.sourceLocal = g.RRis.'*(g.pTx-g.pRis);
g.directionLocal = g.RRis.'*g.dTx;
g.rxDesired = rx.position(:,k); g.rxPlaneZ = rx.nominal(3);
g.hIn = c.ris.thickness*c.optics.screenFraction;
g.hOut = c.ris.thickness-g.hIn; g.zTop = c.ris.thickness;
g.time = m.t(k); g.incidenceDeg = acosd(g.directionLocal(3));
assert(g.sourceLocal(3)<0 && g.directionLocal(3)>0, ...
    'rp:Incidence','TX must be below the lower face and point toward it.');
end
