classdef Fig1PhysicsTest < matlab.unittest.TestCase
    % Independent identities for bare-interface MC and rotated waist geometry.
    methods (Test)
        function rotatedGaussianWaist(t)
            angle=pi/7;R=[cos(angle),-sin(angle);sin(angle),cos(angle)];
            radii=[.4 .1];C=R*diag((radii/2).^2)*R.';
            e=rp_spot_ellipse(struct('covarianceM2',C,'centerXY',[2 -3]));
            t.verifyEqual(e.semiAxesM,radii,'AbsTol',1e-14);
            quadratic=sum((e.offsetXY/C).*e.offsetXY,2);
            t.verifyEqual(quadratic,4*ones(size(quadratic)),'AbsTol',1e-12);
            t.verifyEqual(e.axisEndsXY*e.axisEndsXY.',4*C,'AbsTol',1e-14);
        end
        function normalRayAndFresnelEnergy(t)
            [c,g,mc]=fixture();r=rp_mc_water_air(c,g,mc);
            expected=1-((c.optics.nWater-c.optics.nAir)/(c.optics.nWater+c.optics.nAir))^2;
            t.verifyEqual(r.landingXY,[.1 -.2],'AbsTol',1e-14);
            t.verifyEqual(r.receiverPowerW,2*expected,'AbsTol',1e-14);
            t.verifyEqual(r.receiverPowerW+r.reflectedPowerW,2,'AbsTol',1e-14);
        end
        function obliqueSnellLanding(t)
            [c,g,mc]=fixture();theta=deg2rad(35);
            mc.directionLocal=[sin(theta),0,cos(theta)];r=rp_mc_water_air(c,g,mc);
            airAngle=asin(c.optics.nWater/c.optics.nAir*sin(theta));
            t.verifyEqual(r.landingXY,[.1+g.rxPlaneZ*tan(airAngle),-.2],'AbsTol',1e-13);
            t.verifyEqual(c.optics.nAir*r.directionAir(1,1),c.optics.nWater*sin(theta),'AbsTol',1e-14);
        end
        function totalInternalReflection(t)
            [c,g,mc]=fixture();theta=deg2rad(60);
            mc.directionLocal=[sin(theta),0,cos(theta)];r=rp_mc_water_air(c,g,mc);
            t.verifyEmpty(r.landingXY);t.verifyEqual(r.receiverPowerW,0);
            t.verifyEqual(r.reflectedPowerW,2,'AbsTol',1e-14);
        end
        function noDeviceOrGlassDependence(t)
            [c,g,mc]=fixture();a=rp_mc_water_air(c,g,mc);
            c.optics.nGlass=9;c.array.apertureM=[0 0];c.ris.thickness=999;
            b=rp_mc_water_air(c,g,mc);t.verifyEqual(a,b);
        end
    end
end

function [c,g,mc]=fixture()
c=struct('optics',struct('nWater',1.33,'nAir',1),'source',struct('power',2));
g=struct('pRis',zeros(3,1),'RRis',eye(3),'rxPlaneZ',5);
mc=struct('bottomXY',[.1 -.2],'directionLocal',[0 0 1], ...
    'powerW',2,'ballistic',true,'photons',1);
end
