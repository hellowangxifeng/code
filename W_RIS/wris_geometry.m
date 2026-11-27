function g = wris_geometry(c,m,tx,k,mode)
%WRIS_GEOMETRY Freeze a nanosecond optical snapshot of actual moving bodies.
if nargin<5, mode='joint'; end
g=wris_reference(c);
if ismember(mode,{'joint','ris-only'})
    g.RRis=m.rotationMatrix(:,:,k);
    g.pRis=m.bottomCenter(:,k);
end
if ismember(mode,{'joint','tx-only'})
    g.pTx=tx.position(:,k); g.RTx=tx.rotation(:,:,k); g.dTx=tx.direction(:,k);
end
g.normal=g.RRis(:,3); g.distance=norm(g.pRis-g.pTx);
g.thetaDeg=acosd(min(1,abs(dot(g.dTx,g.normal))));
g.frontSide=dot(g.pRis-g.pTx,g.normal)>0;
[g.axisLocal,g.axisRange,g.axisInside]=wris_intersect(g.pTx.',g.dTx.',g.pRis,g.RRis,[c.ris.Lx c.ris.Ly]);
g.time=m.t(k); g.mode=mode;
% A circular beam's projected major axis lies along its tangential direction.
% This is geometric orientation before aperture clipping, modulo 180 degrees.
localDirection=g.RRis.'*g.dTx;
g.directionLocal=localDirection;
g.yawRisDeg=atan2d(g.RRis(2,1),g.RRis(1,1));
g.projectedOrientationDeg=NaN;
if norm(localDirection(1:2))>1e-10
    g.projectedOrientationDeg=mod(atan2d(localDirection(2),localDirection(1))+90,180)-90;
end
% Remove ONLY yaw as a counterfactual diagnostic, retaining roll/pitch and TX.
noYaw=RIS_rotation([0;0;-deg2rad(g.yawRisDeg)])*g.RRis;
directionNoYaw=noYaw.'*g.dTx;
g.orientationWithoutRisYawDeg=NaN;
if norm(directionNoYaw(1:2))>1e-10
    g.orientationWithoutRisYawDeg=mod(atan2d(directionNoYaw(2),directionNoYaw(1))+90,180)-90;
end
end
