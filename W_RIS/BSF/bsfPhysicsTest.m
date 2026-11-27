classdef bsfPhysicsTest < matlab.unittest.TestCase
    %BSFPHYSICSTEST Physics checks independent of the transport implementation.
    methods (Test)
        function vacuumGaussianFlux(testCase)
            [c,g] = fixture();
            [X,Y] = ndgrid([-0.1,0,0.1],[-0.05,0,0.05]);
            measured = bsf_ballistic(c,g,X,Y);
            width = c.source.divergence*10;
            expected = 2/(pi*width^2)*exp(-2*(X.^2+Y.^2)/width^2);
            testCase.verifyEqual(measured,expected,'RelTol',2e-14);
        end
        function beerLambertUsesRayLength(testCase)
            [c,g] = fixture();
            c.water.absorption = 0.2;
            c.water.scattering = 0.1;
            actual = bsf_ballistic(c,g,0.1,0.2);
            width = c.source.divergence*10;
            expected = 2/(pi*width^2)*exp(-2*0.05/width^2)*exp(-0.3*sqrt(100.05));
            testCase.verifyEqual(actual,expected,'RelTol',2e-14);
        end
        function integratedVacuumPower(testCase)
            [c,g] = fixture();
            x = -0.495:0.01:0.495;
            I = bsf_cell_average(@(X,Y) bsf_ballistic(c,g,X,Y),x,x,6);
            testCase.verifyEqual(sum(I,'all')*0.01^2,1,'AbsTol',1e-12);
        end
        function rigidCoordinateInvariance(testCase)
            [c,g] = fixture();
            expected = bsf_ballistic(c,g,0.07,-0.03);
            angle = 0.43;
            R = [cos(angle),0,sin(angle);0,1,0;-sin(angle),0,cos(angle)];
            shift = [2;-1;3];
            g.pTx = R*g.pTx+shift;
            g.pRis = R*g.pRis+shift;
            g.RTx = R*g.RTx;
            g.RRis = R*g.RRis;
            g.normal = R*g.normal;
            g.dTx = R*g.dTx;
            actual = bsf_ballistic(c,g,0.07,-0.03);
            testCase.verifyEqual(actual,expected,'RelTol',1e-12);
        end
        function noScatteringReturnsZero(testCase)
            [c,g] = fixture();
            out = bsf_scatter_transport(c,g,Photons=1000,Seed=5);
            testCase.verifyEqual(out.Is,zeros(100));
            testCase.verifyEqual(out.power,0);
        end
        function phaseMaskPreservesBallisticPower(testCase)
            [c,g] = fixture();
            c.source.waist = 1.3e-5;
            output = bsf_phase_input(c,g,0.04,0.03,1.7);
            testCase.verifyEqual(abs(output.UafterRis)^2,output.Ib,'RelTol',1e-14);
        end
        function finiteRisAperture(testCase)
            [c,g] = fixture();
            c.source.waist = 1.3e-5;
            output = bsf_phase_input(c,g,0.6,0,0);
            testCase.verifyEqual(output.UafterRis,0);
            testCase.verifyFalse(output.aperture);
        end
    end
end

function [c,g] = fixture()
c = struct('source',struct('power',1,'divergence',0.01,'waist',0, ...
    'wavelength',532e-9,'mediumRefractiveIndex',1.33), ...
    'water',struct('absorption',0,'scattering',0,'anisotropy',0.924), ...
    'mc',struct('maxPath',200,'maxScatter',1000,'chunkSize',10000), ...
    'grid',struct('bins',[100,100]),'ris',struct('Lx',1,'Ly',1));
g = struct('pTx',[0;0;-10],'pRis',[0;0;0],'RTx',eye(3), ...
    'RRis',eye(3),'normal',[0;0;1],'dTx',[0;0;1]);
end
