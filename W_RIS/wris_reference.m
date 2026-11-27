function g = wris_reference(c)
%WRIS_REFERENCE Underwater lower-face optical centre; world Z points upward.
h = RIS_flotation(c.ris);
g.pRis = [c.ris.center;h.glassBottomZ];
angles=c.geometry.rotationXYDeg;
validateattributes(angles,{'numeric'},{'vector','numel',2,'real','finite'});
% Right-hand active fixed-axis rotations: Ry(ay)*Rx(ax)*[0;0;1].
% RIS_rotation([roll;pitch;yaw]) is Rz*Ry*Rx; yaw=0 gives Ry*Rx.
g.RTx=RIS_rotation(deg2rad([angles(:);0]));
g.dTx=g.RTx(:,3);
assert(g.dTx(3)>1e-8,'wris:Direction','Initial beam must point toward the RIS lower face (+z).');
g.pTx = g.pRis-c.geometry.distance*g.dTx;
g.RRis = eye(3); g.normal = [0;0;1];
% Tx optical-frame z is the laser axis, separate from the vehicle yaw.
g.distance = norm(g.pRis-g.pTx);
g.thetaDeg = acosd(abs(dot(g.dTx,g.normal)));
[g.axisLocal,g.axisRange,g.axisInside] = wris_intersect(g.pTx.',g.dTx.',g.pRis,g.RRis,[c.ris.Lx c.ris.Ly]);
end

