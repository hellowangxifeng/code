function result=verify_joint_fading(c)
%VERIFY_JOINT_FADING Independent inverse-CDF QMC vs deterministic EGG/GG solver.
old=rng; clean=onCleanup(@() rng(old)); rng(1211703,'twister');
q=scramble(sobolset(4,'Skip',2048),'MatousekAffineOwen');
u=max(eps,min(1-eps,net(q,2^19))); rows=[]; kinds={'weak','strong'};
power=[1e-4 3e-4 1e-3 3e-3 1e-2 .1];
for ocean=1:2
    for atmosphere=1:2
        cc=c; cc.study.turbulence=kinds{ocean}; cc.study.atmosphere=kinds{atmosphere};
        t=rp_turbulence(cc,32); p=t.oceanParameters; shape=t.ggShapes;
        water=p(4)*gammaincinv(u(:,2),p(3)).^(1/p(5));
        branch=u(:,1)<p(1); water(branch)=-p(2)*log1p(-u(branch,2)); water=water/p(6);
        air=gammaincinv(u(:,3),shape(1))/shape(1).*gammaincinv(u(:,4),shape(2))/shape(2);
        H=water.*air;
        h=struct('gain',1,'fovOutageExact',0,'fovOutageSample',0);
        low=rp_link_metrics(cc,h,power,32); high=rp_link_metrics(cc,h,power,64);
        for k=1:numel(power)
            received=power(k)*H;
            ber=mean(.5*erfc(c.study.responsivityAW*received/(2*sqrt(2)*c.study.noiseStdA)));
            op=mean(received<sqrt(c.study.snrThreshold)*c.study.noiseStdA/c.study.responsivityAW);
            rows=[rows;ocean,atmosphere,power(k),mean(H),low.aber(k),high.aber(k),ber, ...
                low.op(k),high.op(k),op,abs(high.aber(k)-ber),abs(high.op(k)-op)]; %#ok<AGROW>
        end
    end
end
result=array2table(rows,'VariableNames',{'oceanId','atmosphereId','deterministicPowerW','qmcMean', ...
    'berOrder32','berOrder64','berIndependentQMC','opOrder32','opOrder64','opIndependentQMC', ...
    'ber64QMCAbsoluteError','op64QMCAbsoluteError'});
writetable(result,fullfile(c.root,'audit','stage2','joint_fading_verification.csv'));
save(fullfile(c.root,'audit','stage2','joint_fading_verification.mat'),'result');
disp(result);
end
