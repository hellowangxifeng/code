classdef RisLinkPhysicsTest < matlab.unittest.TestCase
    % Independent identities for finite-aperture channel and coherent pixels.
    methods (Test)
        function quantizedCommandLimits(t)
            c=rp_study_config();
            [b,angle]=rp_command(c,[10 -10;0.1 -.1]);
            limit=floor(c.device.maxAngleDeg/c.device.angleStepDeg)*c.device.angleStepDeg;
            t.verifyLessThanOrEqual(max(abs(angle),[],'all'),limit);
            t.verifyEqual(angle/c.device.angleStepDeg,round(angle/c.device.angleStepDeg),'AbsTol',1e-10);
            t.verifyEqual(b,c.optics.nAir*sind(angle),'AbsTol',1e-15);
        end
        function zeroSignal(t)
            c=rp_study_config();
            h=struct('gain',zeros(32,1),'fovOutageExact',1,'fovOutageSample',1);
            result=rp_link_metrics(c,h,[.01 .1 1]);
            t.verifyEqual(result.aber,[.5 .5 .5],'AbsTol',1e-12);
            t.verifyEqual(result.op,[1 1 1],'AbsTol',1e-12);
            t.verifyEqual(result.meanReceivedPowerW,[0 0 0]);
        end
        function deterministicOOK(t)
            c=rp_study_config(); c.study.turbulence='none'; c.study.atmosphere='none';
            h=struct('gain',.01*ones(16,1),'fovOutageExact',0,'fovOutageSample',0);
            m=rp_link_metrics(c,h,.2);
            expected=.5*erfc(c.study.responsivityAW*.002/(2*sqrt(2)*c.study.noiseStdA));
            t.verifyEqual(m.aber,expected,'AbsTol',1e-14);
            t.verifyEqual(m.op,0);
        end
        function fullRayleighTails(t)
            c=rp_study_config(); s=rp_link_samples(c,32768,752);
            r=vecnorm(s.positionNormal,2,2);
            t.verifyEqual(mean(r.^2),2,'AbsTol',.003);
            t.verifyEqual(mean(r>3),exp(-4.5),'AbsTol',.001);
        end
        function gaussianApertureIntegral(t)
            c=rp_study_config(); c.study.positionSigmaM=0; c.study.angleSigmaRad=0;
            c.study.angleBiasRad=[0 0];
            axis=linspace(-.3,.3,1001); [X,Y]=meshgrid(axis);
            width=.1; I=2/(pi*width^2)*exp(-2*(X.^2+Y.^2)/width^2);
            w=struct('x',axis,'y',axis,'irradianceWm2',I);
            g=struct('pRis',[0;0;-.01],'rxPlaneZ',5);
            h=rp_aperture_channel(c,g,w,[0 0],rp_link_samples(c,16));
            expected=1-exp(-2*c.rx.apertureRadius^2/width^2);
            t.verifyEqual(h.gain,expected*ones(16,1),'AbsTol',2e-5);
        end
        function tiltedUniformDisk(t)
            c=rp_study_config(); c.study.positionSigmaM=0; c.study.angleSigmaRad=0;
            c.study.angleBiasRad=[.03 0]; c.study.fovHalfAngleRad=.1;
            axis=linspace(-5,5,101);
            w=struct('x',axis,'y',axis,'irradianceWm2',2*ones(101));
            g=struct('pRis',[0;0;-.01],'rxPlaneZ',5);
            h=rp_aperture_channel(c,g,w,[3 0],rp_link_samples(c,16));
            expected=2*pi*c.rx.apertureRadius^2*cos(.03)/h.nominalDirection(3);
            t.verifyEqual(h.gain,expected*ones(16,1),'AbsTol',1e-12);
        end
        function rayleighFOVFloor(t)
            c=rp_study_config(); c.study.angleBiasRad=[0 0]; c.study.fovHalfAngleRad=.07;
            axis=linspace(-1,1,41);
            w=struct('x',axis,'y',axis,'irradianceWm2',ones(41));
            g=struct('pRis',[0;0;-.01],'rxPlaneZ',5);
            h=rp_aperture_channel(c,g,w,[0 0],rp_link_samples(c,32768,753));
            expected=exp(-c.study.fovHalfAngleRad^2/(2*c.study.angleSigmaRad^2));
            t.verifyEqual(h.fovOutageExact,expected,'AbsTol',1e-12);
            t.verifyEqual(h.fovOutageSample,expected,'AbsTol',.001);
        end
        function pixelInterference(t)
            fields=[1 1i]; slopesX=[0 200]; slopesY=[0 -100]; pixel=[.01 .02];
            value=rp_pixel_intensity(fields,slopesX,slopesY,pixel);
            t.verifyEqual(value,2,'AbsTol',1e-14);
            valueInPhase=rp_pixel_intensity([1 1],slopesX,slopesY,pixel);
            expected=2+2*(sin(1)/1)^2;
            t.verifyEqual(valueInPhase,expected,'AbsTol',1e-12);
        end
        function layeredNormalPath(t)
            path=rp_layer_path([0 0],[.005 5],[1.5 1]);
            t.verifyEqual(path.q,[0 0]);
            t.verifyEqual(path.opl,5.0075,'AbsTol',1e-12);
        end
        function powerAndNoiseMonotonicity(t)
            c=rp_study_config();
            h=struct('gain',[0;.001;.003;.01],'fovOutageExact',.25,'fovOutageSample',.25);
            m=rp_link_metrics(c,h,[.01 .1 1]);
            t.verifyLessThanOrEqual(diff(m.aber),[0 0]);
            t.verifyLessThanOrEqual(diff(m.op),[0 0]);
            t.verifyGreaterThanOrEqual(m.aber,.125*ones(1,3));
            t.verifyGreaterThanOrEqual(m.op,.25*ones(1,3));
        end
    end
end
