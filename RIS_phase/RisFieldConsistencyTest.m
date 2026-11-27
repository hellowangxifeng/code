classdef RisFieldConsistencyTest < matlab.unittest.TestCase
    % Coordinate, hybrid-scatter and electric/flux normalization identities.
    properties (TestParameter)
        pose = struct('translated',eye(3),'tilted', ...
            [cos(.3) 0 sin(.3);0 1 0;-sin(.3) 0 cos(.3)])
        transverseQ = struct('normal',0,'oblique',.6,'nearCritical',.97)
    end
    methods (Test)
        function incidentRigidPose(t,pose)
            c=rp_study_config();
            base=localScene([0;0;0],eye(3));
            moved=localScene([12;-7;3],pose);
            xy=[0 0;.001 -.002;-.003 .001];
            first=rp_incident_field(c,base,xy);
            second=rp_incident_field(c,moved,xy);
            t.verifyEqual(moved.sourceLocal,base.sourceLocal,'AbsTol',2e-14);
            t.verifyEqual(second.q,first.q,'AbsTol',2e-14);
            t.verifyEqual(second.irradianceWm2,first.irradianceWm2,'RelTol',1e-11);
            t.verifyEqual(second.field,first.field,'RelTol',2e-7);
        end
        function receiverTranslation(t)
            c=rp_study_config(); shift=[12;-7;3];
            base=localScene([0;0;0],eye(3));
            moved=localScene(shift,eye(3));
            reference=[.01 -.02]; receiver=[2 -1;2.1 -.9];
            first=rp_paraxial(c,base,reference,receiver);
            second=rp_paraxial(c,moved,reference,receiver+shift(1:2).');
            t.verifyEqual(second.q,first.q,'AbsTol',2e-14);
            t.verifyEqual(second.opl,first.opl,'AbsTol',2e-13);
            t.verifyEqual(second.B,first.B,'AbsTol',2e-13);
            t.verifyTrue(all(second.valid));
        end
        function scatteredBottomToLC(t)
            [c,g,mc]=scatterInput();
            a=rp_illuminate(c,g,mc);
            expectedQ=c.optics.nWater*mc.directionLocal(:,1:2);
            expectedPower=mc.powerW.*rp_fresnel(c.optics.nWater,c.optics.nGlass,expectedQ);
            t.verifyEqual(a.lcXY,zeros(3,2),'AbsTol',1e-15);
            t.verifyEqual(a.qIn,expectedQ,'AbsTol',1e-15);
            t.verifyEqual(a.powerW,expectedPower,'AbsTol',1e-15);
            t.verifyFalse(any(a.ballistic));
            t.verifyEqual(a.activeDirectPowerW,0,'AbsTol',1e-15);
        end
        function scatteredScreenToReceiver(t)
            [c,g,mc]=scatterInput(); a=rp_illuminate(c,g,mc);
            p=scatterPhase(c); r=rp_trace(c,g,a,p);
            expectedQ=a.qIn+p.gradientOverK0;
            topXY=g.hOut*expectedQ(1:2,:)./sqrt(c.optics.nGlass^2-sum(expectedQ(1:2,:).^2,2));
            expectedDirection=[expectedQ(1:2,:),sqrt(1-sum(expectedQ(1:2,:).^2,2))];
            expectedLanding=[topXY,repmat(g.zTop,2,1)]+ ...
                (g.rxPlaneZ-g.zTop)*expectedDirection./expectedDirection(:,3);
            expectedPower=a.powerW.*rp_fresnel(c.optics.nGlass,c.optics.nAir,expectedQ);
            t.verifyEqual(r.qOut,expectedQ,'AbsTol',1e-15);
            t.verifyEqual(r.valid,[true;true;false]);
            t.verifyEqual(r.landing(1:2,:),expectedLanding,'AbsTol',2e-13);
            t.verifyTrue(all(isnan(r.landing(3,:))));
            t.verifyEqual(r.powerW,expectedPower,'AbsTol',1e-15);
            t.verifyEqual(r.powerW(3),0,'AbsTol',1e-15);
        end
        function incoherentScatterPixels(t)
            [c,g,mc]=scatterInput(); a=rp_illuminate(c,g,mc);
            p=scatterPhase(c); r=rp_trace(c,g,a,p);
            w=rp_wave(c,g,a,mc,p,r);
            expectedPower=sum(r.powerW(r.valid));
            [~,ix]=min(abs(w.x-r.landing(1,1)));
            [~,iy]=min(abs(w.y-r.landing(1,2)));
            t.verifyEqual(w.scatteredMapPowerW,expectedPower,'AbsTol',1e-14);
            t.verifyEqual(w.scatteredIrradianceWm2(iy,ix)*w.cellAreaM2,expectedPower,'AbsTol',1e-14);
            t.verifyEqual(w.irradianceWm2,w.scatteredIrradianceWm2,'AbsTol',1e-14);
            t.verifyEqual(w.ballisticPixelMeanWm2,zeros(size(w.irradianceWm2)),'AbsTol',1e-14);
        end
        function electricAndFluxTransmission(t,transverseQ)
            c=rp_study_config(); ng=c.optics.nGlass; na=c.optics.nAir;
            kg=sqrt(ng^2-transverseQ^2); ka=sqrt(na^2-transverseQ^2);
            epsilon0=8.8541878128e-12; lightSpeed=299792458;
            incidentElectric=2; factor=epsilon0*lightSpeed/2;
            ts=2*kg/(kg+ka);
            tp=2*ng*na*kg/(na^2*kg+ng^2*ka);
            inputFlux=factor*kg*incidentElectric^2;
            outputFlux=factor*ka*incidentElectric^2*(ts^2+tp^2)/2;
            transmission=rp_fresnel(ng,na,[transverseQ 0]);
            t.verifyEqual(outputFlux,inputFlux*transmission,'AbsTol',1e-15);
            t.verifyEqual(sqrt(outputFlux),sqrt(inputFlux)*sqrt(transmission),'AbsTol',1e-14);
        end
    end
end

function g=localScene(origin,rotation)
source=[-2;1;-10]; direction=-source/norm(source);
first=cross([0;1;0],direction); first=first/norm(first);
basis=[first,cross(direction,first),direction];
m=struct('rotationMatrix',rotation,'bottomCenter',origin,'t',0);
tx=struct('position',origin+rotation*source,'direction',rotation*direction, ...
    'rotation',rotation*basis);
rx=struct('position',origin+[0;0;5.01],'nominal',origin+[0;0;5.01]);
g=rp_geometry(rp_study_config(),m,tx,rx,1);
end

function [c,g,mc]=scatterInput()
c=rp_study_config(); c.array.count=[1 1]; c.source.power=0;
c.spot.mapBins=[16 16];
g=localScene([0;0;0],eye(3));
q=[.2 -.1;.2 -.1;.9 .1];
bottom=-g.hIn*q./sqrt(c.optics.nGlass^2-sum(q.^2,2));
direction=[q,sqrt(c.optics.nWater^2-sum(q.^2,2))]/c.optics.nWater;
mc=struct('bottomXY',bottom,'directionLocal',direction,'powerW',[.1;.2;.3], ...
    'ballistic',false(3,1));
end

function p=scatterPhase(c)
gradient=rp_command(c,[.15 0]);
p=struct('referenceXY',[0 0],'gradientOverK0',gradient,'phaseCurvature',zeros(2), ...
    'phaseModel','linear','phaseAtReferenceRad',0,'unitIds',1,'centerXY',[2 -1], ...
    'receiverHalfWidthM',[3 3],'outgoing',struct('valid',true));
end
