classdef WrdPhysicsTest < matlab.unittest.TestCase
    %WRDPHYSICSTEST Physics/bookkeeping regressions for RIS_dynamic.
    %   results = runtests('RIS_dynamic/tests/WrdPhysicsTest.m')
    %   None of these tests solves the six-DOF RIS motion; together they take
    %   well under a minute.

    methods (TestClassSetup)
        function addPaths(~)
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(fileparts(here),'algorithm'));
            wrd_setup();
        end
    end

    methods (Test)
        function matchesFrozenFourierField(testCase)
            % Interpolated height field and slopes equal the exact Fourier sum.
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

        function flatSurfaceIsTheStillWaterPlane(testCase)
            c = wrd_config('smoke');
            sea = wrd_flat_surface();
            testCase.verifyEqual(sea.amplitudeBound,0);
            [eta,gradient] = sea.evaluate([0 0; 3 -4]);
            testCase.verifyEqual(eta,[0;0]);
            testCase.verifyEqual(gradient,zeros(2,2));
            distance = wrd_crossing([0.2 0.1 -1],[0 0 1],5,sea,c.surface,1);
            testCase.verifyEqual(distance,1,'AbsTol',1e-12);
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
            % A 15-degree facet turns a 55-degree (TIR) ray into a 40-degree one.
            [c,g,sea] = collimatedCase(55);
            sea.amplitudeBound = 10;
            sea.slopeBound = [tand(15) 0];
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

        function staticRisLinkReachesFixedReceiver(testCase)
            % Still RIS at its static draft: energy closes, the steered beam
            % reaches the common RX and nothing leaves under the flat sea.
            c = wrd_config('smoke');
            c.geometry.rotationXYDeg = [0 40];
            ref = wris_reference(c);
            m = struct('rotationMatrix',eye(3),'bottomCenter',ref.pRis,'t',0);
            tx = struct('position',ref.pTx,'direction',ref.dTx,'rotation',ref.RTx);
            rx = struct('position',c.receiver.center(:),'nominal',c.receiver.center(:));
            g = rp_geometry(c,m,tx,rx,1);
            out = wrd_ris_mc(c,g,c.mc.seed,wrd_flat_surface());
            testCase.verifyLessThan(abs(out.energyResidual),1e-9);
            testCase.verifyGreaterThan(out.rxPowerW,0);
            testCase.verifyLessThanOrEqual(out.rxPowerW,out.transmittedPowerW*(1+1e-12));
            testCase.verifyLessThanOrEqual(out.transmittedPowerW,out.incidentPowerW*(1+1e-12));
            testCase.verifyEqual(out.overwashedExitPowerW,0);
        end

        function summaryStatistics(testCase)
            c = wrd_config('smoke');
            c.metrics.outageThresholdW = 1e-6;
            direct = [fakeLink(0,0),fakeLink(2e-6,5),fakeLink(4e-6,5)];
            ris = [fakeLink(1e-3,9),fakeLink(1e-3,9),fakeLink(1e-3,9)];
            frames = struct('t',{0,0.5,1},'absoluteTime',{25,25.5,26},'etaM',{0,0,0}, ...
                'slope',{[0 0],[0 0],[0 0]},'risTiltDeg',{0,0,0},'bottomDryFraction',{0,0,0}, ...
                'topWetFraction',{0,0,0},'direct',num2cell(direct),'flat',num2cell(direct), ...
                'ris',num2cell(ris));
            item = struct('windMps',3,'windDirectionDeg',50,'txRotationXYDeg',[0 40], ...
                'nominalIncidenceDeg',40,'significantWaveHeightM',0.1,'frames',frames(:));
            [frameTable,caseTable] = wrd_summary({item},c);
            testCase.verifyEqual(height(frameTable),3);
            rough = caseTable(caseTable.link=="roughDirect",:);
            testCase.verifyEqual(rough.meanW,2e-6,'AbsTol',1e-18);
            testCase.verifyEqual(rough.outageProbability,1/3,'AbsTol',1e-12);
            testCase.verifyEqual(rough.zeroHitFraction,1/3,'AbsTol',1e-12);
            testCase.verifyEqual(rough.scintillationIndex,2/3,'AbsTol',1e-12);
            testCase.verifyEqual(rough.risOverRoughDirectDb,10*log10(500),'AbsTol',1e-10);
            testCase.verifyEqual(height(caseTable),3);
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
sea = wrd_flat_surface();
end

function [eta,gradient] = slopingSurface(xy)
gradient = repmat([-tand(15),0],size(xy,1),1);
eta = xy(:,1).*gradient(:,1);
end

function link = fakeLink(power,hits)
link = struct('rxPowerW',power,'powerSEW',0,'hitCount',hits,'ballisticPowerW',power, ...
    'scatterPowerW',0,'rxOnlyApertureW',power,'tirPowerW',0,'incidentPowerW',1, ...
    'transmittedPowerW',1,'waveBlockedPowerW',0,'energyResidual',0,'activePowerW',1, ...
    'saturatedPowerFraction',0,'overwashedExitPowerW',0);
end
