classdef RdPhysicsTest < matlab.unittest.TestCase
    %RDPHYSICSTEST Fast physics / bookkeeping checks of RIS_dynamic (no RIS_float, no wave run).
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
            air = rd_direct_mc(c,g,flatSurface(),7);
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
            air = rd_direct_mc(c,g,flatSurface(),7);
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

        function exactFovTailInMetrics(tc)
            % Samples outside the FOV are replaced by the exact Rice mass p_F.
            c = tc.c;
            G = 0.01*rand(64,1);
            psi = [zeros(60,1);ones(4,1)];
            P = [1e-3 1e-1 10];
            m = rd_link_metrics(c,G,psi,P,'weak','weak');
            cc = c; cc.study.turbulence = 'weak'; cc.study.atmosphere = 'weak';
            r = rp_link_metrics(cc,struct('gain',G(1:60),'fovOutageExact',0,'fovOutageSample',0), ...
                P,c.metrics.order);
            pF = m.fovOutageExact;
            tc.verifyGreaterThan(pF,0);
            tc.verifyEqual(m.aber,(1-pF)*r.aber+pF/2,'AbsTol',1e-15);
            tc.verifyEqual(m.op,(1-pF)*r.op+pF,'AbsTol',1e-15);
        end

        function neighbourCoherenceMatchesRisPhase(tc)
            % With every pair inside the neighbour rings, the result equals
            % RIS_phase/rp_pixel_intensity on the same dense fields.
            rng(3);
            pixels = 40; K = 9;
            fields = complex(randn(pixels,K),randn(pixels,K));
            sx = 400*randn(pixels,K); sy = 400*randn(pixels,K);
            [i,j] = ndgrid(1:3,1:3);
            support = repmat({(1:pixels).'},1,K);
            values = num2cell(fields,1); slopeX = num2cell(sx,1); slopeY = num2cell(sy,1);
            I = rd_pixel_intensity(support,values,slopeX,slopeY,pixels,[0.01 0.01],[i(:),j(:)],2);
            I0 = rp_pixel_intensity(fields,sx,sy,[0.01 0.01],0);
            tc.verifyEqual(I,I0,'RelTol',1e-10,'AbsTol',1e-12);
        end

        function selectionRules(tc)
            % Flat and Robust minimise their own family; Max mean power maximises <P_r>.
            family = ["steering","flat_top","flat_top","robust","robust"];
            A = [.4 .4; .1 .2; .2 .05; .3 .01; .05 .3];
            O = A; M = repmat([1;2;3;4;5],1,2);
            ch = rd_select(A,O,M,family,true(1,5),'per_power',1);
            tc.verifyEqual(squeeze(ch(2,:,1)),[1 1]);
            tc.verifyEqual(squeeze(ch(3,:,1)),[2 3]);
            tc.verifyEqual(squeeze(ch(4,:,1)),[5 4]);
            tc.verifyEqual(squeeze(ch(5,:,1)),[5 5]);
            frozen = rd_select(A,O,M,family,true(1,5),'design_power',2);
            tc.verifyEqual(squeeze(frozen(3,:,1)),[3 3]);
        end

        function requiredPowerRanking(tc)
            c = tc.c;
            P = 0:1:40;
            base = @(shift) min(0.5,10.^(-(P-shift)/5));
            curves.aber = [0.5*ones(1,numel(P));base(20);base(12);base(5);base(10)];
            curves.op = curves.aber;
            r = rd_rank(curves,P,c);
            tc.verifyEqual(r.aber.winners,4);
            tc.verifyTrue(r.aber.clear);
            tc.verifyEqual(r.aber.requiredDBm(1),Inf);
            tc.verifyEqual(r.aber.gainOverSteeringDB,15,'AbsTol',1e-9);
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

function sea = flatSurface()
% Still-water plane z = 0 through the same rough-surface interface code.
motion.parameters.startTime = 0;
motion.wave = struct('omega',0,'coefficient',0,'kx',0,'ky',0);
settings = struct('gridStep',1,'gridHalfWidth',0,'marchStep',1,'rootTolerance',1e-9);
sea = rd_surface(motion,0,settings,[0;0]);
end

function m = toySea()
m.parameters.startTime = 0;
m.wave = struct('omega',[1.1;2.0;3.1],'kx',[0.6;1.5;-2.4],'ky',[0.2;-0.9;1.7], ...
    'coefficient',[0.05;0.02i;0.01-0.005i]);
end
