classdef RdPhysicsTest < matlab.unittest.TestCase
    %RDPHYSICSTEST Fast physics / bookkeeping checks of RIS_dynamic (no RIS_float).
    %   r = runtests('tests/RdPhysicsTest.m'); table(r)
    properties
        c
    end

    methods (TestClassSetup)
        function setup(tc)
            here = fileparts(fileparts(mfilename('fullpath')));
            addpath(fullfile(here,'algorithm'),fullfile(here,'plots'));
            rd_setup();
            tc.c = rd_config('smoke');
        end
    end

    methods (Test)
        function parallelBeamCapture(tc)
            % Uniform horizontal flux I_h along d, disc normal = d:
            % captured = I_h * pi a^2 / d_z (RIS_phase projected-area rule).
            c = tc.c;
            theta = deg2rad(30);
            d = [sin(theta) 0 cos(theta)];
            [pk,Ih] = uniformBeam(d,1e-3,0.15);
            rx = rd_receiver(c,[0 0 5]-5/d(3)*d,[0 0],5,zeroSamples(4),pointingFree(c,0.1));
            gain = rd_capture(c,pk,rx);
            expected = Ih*pi*c.rx.apertureRadius^2/d(3)/c.source.power;
            tc.verifyEqual(gain,repmat(expected,4,1),'RelTol',0.01);
        end

        function perRayFovRejectsOffAxisLight(tc)
            c = tc.c;
            d = [sin(0.2) 0 cos(0.2)];
            pk = uniformBeam(d,2e-3,0.2);
            rx = rd_receiver(c,[0 0 0],[0 0],5,zeroSamples(2),pointingFree(c,0.1));
            tc.verifyEqual(rd_capture(c,pk,rx),zeros(2,1));
            rx = rd_receiver(c,[0 0 0],[0 0],5,zeroSamples(2),pointingFree(c,0.25));
            tc.verifyGreaterThan(min(rd_capture(c,pk,rx)),0);
        end

        function jitterMovesTheDisc(tc)
            % An offset of 3.5a puts the disc completely outside a 4a-wide strip.
            c = tc.c;
            a = c.rx.apertureRadius;
            [pk,~] = uniformBeam([0 0 1],1e-3,0.12);
            keep = abs(pk.xy(:,1))<=2*a;
            pk = struct('xy',pk.xy(keep,:),'dir',pk.dir(keep,:),'power',pk.power(keep), ...
                'ballistic',pk.ballistic(keep));
            s = zeroSamples(2);
            s.positionNormal = [0 0; 1 0];
            o = pointingFree(c,0.1);
            o.positionSigmaM = 3.5*a;
            rx = rd_receiver(c,[0 0 0],[0 0],5,s,o);
            gain = rd_capture(c,pk,rx);
            tc.verifyGreaterThan(gain(1),0);
            tc.verifyEqual(gain(2),0);
        end

        function flatSurfaceSnellFresnel(tc)
            c = clearWater(tc.c);
            g = sourceGeometry(40);
            air = rd_direct_mc(c,g,rd_flat_surface(),7);
            q = c.optics.nWater*sind(40);
            tc.verifyEqual(vecnorm(air.dir(:,1:2),2,2),repmat(q/c.optics.nAir,numel(air.power),1), ...
                'AbsTol',1e-4);
            T = rp_fresnel(c.optics.nWater,c.optics.nAir,[q 0]);
            tc.verifyEqual(sum(air.power),T*c.source.power,'RelTol',1e-4);
            tc.verifyLessThan(abs(air.energyResidual),1e-12);
        end

        function flatSurfaceTotalInternalReflection(tc)
            c = clearWater(tc.c);
            g = sourceGeometry(55);
            air = rd_direct_mc(c,g,rd_flat_surface(),7);
            tc.verifyEmpty(air.power);
            tc.verifyEqual(air.ballisticTirFraction,1,'AbsTol',1e-12);
            tc.verifyLessThan(abs(air.energyResidual),1e-12);
        end

        function roughSurfaceMatchesWaveField(tc)
            c = tc.c;
            m = toySea();
            sea = rd_surface(m,1.3,c.surface,[0 0]);
            xy = [0.4*randn(200,2);[3 -2.5;-4 1]];
            [eta,slope] = sea.evaluate(xy);
            [eta0,slope0] = RIS_wave_field(m,xy,1.3);
            tc.verifyEqual(eta,eta0,'AbsTol',1e-5);
            tc.verifyEqual(slope,slope0,'AbsTol',1e-3);
        end

        function roughSurfaceLedgerCloses(tc)
            c = tc.c;
            c.mc.photons = 2000; c.mc.chunkSize = 1000;
            g = sourceGeometry(40);
            sea = rd_surface(toySea(),0.7,c.surface,[0 0]);
            air = rd_direct_mc(c,g,sea,11);
            tc.verifyLessThan(abs(air.energyResidual),1e-9);
            tc.verifyLessThanOrEqual(sum(air.power),air.transmittedPowerW+1e-12);
        end

        function compressionConservesPowerAndCentroid(tc)
            n = 5000;
            d = [0.002*randn(n,2),ones(n,1)];
            pk = struct('xy',0.005*randn(n,2),'dir',d./vecnorm(d,2,2),'power',rand(n,1), ...
                'ballistic',rand(n,1)>0.3);
            out = rd_compress(pk,struct('positionM',5e-3,'directionCos',5e-3));
            tc.verifyLessThan(numel(out.power),n);
            tc.verifyEqual(sum(out.power),sum(pk.power),'RelTol',1e-12);
            tc.verifyEqual(out.power.'*out.xy,pk.power.'*pk.xy,'AbsTol',1e-12);
            tc.verifyEqual(vecnorm(out.dir,2,2),ones(numel(out.power),1),'AbsTol',1e-12);
        end

        function timeAverageOfFrameMetrics(tc)
            c = tc.c;
            G = 0.01*rand(64,3);
            G(1:5,2) = 0;
            rx = struct('fovOutageExact',0,'fovOutageSample',0);
            P = [1e-3 1e-2 1];
            m = rd_metrics(c,G,rx,P);
            aber = zeros(3,numel(P));
            for f = 1:3
                r = rp_link_metrics(c,struct('gain',G(:,f),'fovOutageExact',0,'fovOutageSample',0),P);
                aber(f,:) = r.aber;
            end
            tc.verifyEqual(m.aber,mean(aber,1),'AbsTol',1e-15);
            tc.verifyEqual(m.meanGain,mean(G(:)),'RelTol',1e-12);
        end
    end
end

function [pk,Ih] = uniformBeam(d,h,half)
[x,y] = meshgrid(-half:h:half);
n = numel(x);
Ih = 1;
pk = struct('xy',[x(:),y(:)],'dir',repmat(d/norm(d),n,1),'power',repmat(Ih*h^2,n,1), ...
    'ballistic',true(n,1));
end

function s = zeroSamples(n)
s = struct('positionNormal',zeros(n,2),'angleNormal',zeros(n,2),'n',n,'seed',0);
end

function o = pointingFree(c,fov)
o = c.study;
o.positionSigmaM = 0;
o.angleSigmaRad = 0;
o.angleBiasRad = [0 0];
o.fovHalfAngleRad = fov;
end

function c = clearWater(c)
c.water.absorption = 0;
c.water.scattering = 0;
c.mc.photons = 200;
c.mc.chunkSize = 200;
c.source.divergence = 1e-7;
c.source.waist = 1e-4;
end

function g = sourceGeometry(incidenceDeg)
R = [cosd(incidenceDeg) 0 sind(incidenceDeg); 0 1 0; -sind(incidenceDeg) 0 cosd(incidenceDeg)];
g = struct('RTx',R,'dTx',R(:,3),'pTx',[0;0;-5]);
end

function m = toySea()
m.parameters.startTime = 0;
m.wave = struct('omega',[1.1;2.0;3.1],'kx',[0.6;1.5;-2.4],'ky',[0.2;-0.9;1.7], ...
    'coefficient',[0.05;0.02i;0.01-0.005i]);
end
