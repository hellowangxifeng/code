function result=verify_direct_incidence(c)
%VERIFY_DIRECT_INCIDENCE Check TIR tails and independent direct-beam energy.
rows=[];
for angle=c.study.angleGridDeg
    [g,~,~,cc]=rp_stage2_scene(c,7.5,angle);
    for z=[7.5 10]
        g.rxPlaneZ=g.pRis(3)+c.ris.thickness+z;
        b=rp_direct_wave(cc,g); m=b.designMetrics;
        capture=b.spot.powerW/max(b.wave.angularCubaturePowerW,realmin);
        rows=[rows;angle,z,b.wave.angularCubaturePowerW,b.spot.powerW,capture, ...
            m.meanReceivedPowerW,m.aber,m.op]; %#ok<AGROW>
    end
end
result=array2table(rows,'VariableNames',{'incidenceDeg','heightM','angularPowerW', ...
    'mapPowerW','mapToAngularRatio','meanAperturePowerW','aber','op'});
writetable(result,fullfile(c.root,'audit','stage2','direct_incidence_verification.csv'));
save(fullfile(c.root,'audit','stage2','direct_incidence_verification.mat'),'result');
disp(result);
end
