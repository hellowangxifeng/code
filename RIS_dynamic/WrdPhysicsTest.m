classdef WrdPhysicsTest < matlab.unittest.TestCase
    methods (Test)
        function matchesFrozenFourierField(testCase)
            c = wrd_config('smoke');
            m.wave = RIS_wave_components(c.ris);
            m.parameters = c.ris;
            c.surface.gridHalfWidth = 0.3;
            sea = wrd_surface(m,0.3,c.surface,[0;0]);
            xy = [0 0; 0.075 -0.075; 1 -2];
            [expected,slope] = RIS_wave_field(m,xy,0.3);
            [height,gradient] = sea.evaluate(xy);
            testCase.verifyEqual(height,expected,'AbsTol',1e-12);
            testCase.verifyEqual(gradient(3,:),slope(3,:),'AbsTol',1e-12);
        end

        function findsFirstOfMultipleCrossings(testCase)
            c = wrd_config('smoke');
            m.parameters.startTime = 0;
            m.wave = struct('omega',0,'coefficient',1,'kx',1,'ky',0);
            c.surface.gridHalfWidth = 0;
            sea = wrd_surface(m,0,c.surface,[0;0]);
            distance = wrd_crossing([0 0 0],[1 0 0],2*pi,sea,c.surface,1);
            testCase.verifyEqual(distance,pi/2,'AbsTol',3e-6);
            recross = wrd_crossing([pi 0 0],[1 0 0],2*pi,sea,c.surface,-1);
            testCase.verifyEqual(recross,pi/2,'AbsTol',3e-6);
        end

        function flatWaterMatchesFresnelAndSnell(testCase)
            [c,g,sea] = collimatedCase(40);
            ray = [c.optics.nWater*sind(40),0,sqrt(1-(c.optics.nWater*sind(40))^2)];
            c.receiver.center = ray/ray(3);
            c.receiver.normal = ray;
            result = wrd_direct_mc(c,g,sea,1);
            expected = rp_fresnel(c.optics.nWater,c.optics.nAir,[ray(1),0]);
            testCase.verifyEqual(result.rxPowerW,c.source.power*expected,'AbsTol',1e-12);
            testCase.verifyEqual(result.hitCount,c.mc.photons);
            testCase.verifyEqual(result.energyResidual,0,'AbsTol',1e-12);
        end

        function flatWaterAboveCriticalIsTotalReflection(testCase)
            [c,g,sea] = collimatedCase(55);
            result = wrd_direct_mc(c,g,sea,1);
            testCase.verifyEqual(result.tirPowerW,c.source.power,'AbsTol',1e-12);
            testCase.verifyEqual(result.transmittedPowerW,0);
            testCase.verifyEqual(result.energyResidual,0,'AbsTol',1e-12);
        end

        function slopeOpensAboveCriticalRay(testCase)
            [c,g,sea] = collimatedCase(55);
            tilt = deg2rad(15);
            sea.amplitudeBound = 10;
            sea.slopeBound = [tan(tilt) 0];
            sea.evaluate = @slopingSurface;
            expected = rp_fresnel(c.optics.nWater,c.optics.nAir,[c.optics.nWater*sind(40),0]);
            result = wrd_direct_mc(c,g,sea,1);
            testCase.verifyEqual(result.transmittedPowerW,c.source.power*expected,'AbsTol',1e-12);
            testCase.verifyEqual(result.tirPowerW,0);
        end

        function rejectsOutsideFovAndAperture(testCase)
            c = wrd_config('smoke');
            hits = [0 0 1 1 1; 0 0 1 0 cosd(20); 0.1 0 1 0 1];
            result = wrd_score(c,hits);
            scale = c.source.power/c.mc.photons;
            testCase.verifyEqual(result.rxPowerW,scale,'AbsTol',eps);
            testCase.verifyEqual(result.rxOnlyApertureW,2*scale,'AbsTol',eps);
            testCase.verifyEqual(result.mapPowerW,3*scale,'AbsTol',eps);
        end
    end
end

function [c,g,sea] = collimatedCase(angle)
c = wrd_config('smoke');
c.water.absorption = 0;
c.water.scattering = 0;
c.source.waist = 0;
c.source.divergence = 0;
g.RTx = RIS_rotation([0;deg2rad(angle);0]);
g.dTx = g.RTx(:,3);
g.pTx = -2*g.dTx;
m.parameters.startTime = 0;
m.wave = struct('omega',0,'coefficient',0,'kx',0,'ky',0);
sea = wrd_surface(m,0,c.surface,[0;0]);
end

function [eta,gradient] = slopingSurface(xy)
gradient = repmat([-tand(15),0],size(xy,1),1);
eta = xy(:,1).*gradient(:,1);
end
