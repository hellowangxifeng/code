function validation = bsf_validate_outputs()
%BSF_VALIDATE_OUTPUTS Independent seeds, packet accounting, and paper quadrature.
root = bsf_setup();
folder = fullfile(root,'results');
loaded = load(fullfile(fileparts(root),'results','wris_results.mat'),'s');
c = loaded.s.config;
panels = [1,1;3,3;6,3];
records = cell(size(panels,1),1);
for j = 1:size(panels,1)
    row = panels(j,1);
    col = panels(j,2);
    saved = load(fullfile(folder,sprintf('panel_r%d_c%d.mat',row,col)),'result');
    r = saved.result;
    reference = loaded.s.panels{row,col}.out;
    independent = bsf_scatter_transport(c,r.geometry,Photons=8e6,Seed=710000+j);
    difference = independent.Is-r.Is;
    predictedSquaredNoise = sum(independent.standardError.^2+r.standardError.^2,'all');
    ratio = norm(difference,'fro')/sqrt(predictedSquaredNoise);
    powerZ = (independent.power-r.powerMean)/hypot(independent.powerSE,r.powerSE);
    records{j} = struct('row',row,'column',col,'independentPhotons',independent.photons, ...
        'errorOverPredictedRMS',ratio,'independentPowerZ',powerZ, ...
        'compatible',abs(powerZ)<5 && ratio>0.8 && ratio<1.2);
    fprintf('Independent r%d c%d: RMS/SE=%.5f, power z=%.4f\n',row,col,ratio,powerZ);
    if j==1
        [replay,packets] = wris_mc(c,r.geometry,reference.seed,1);
        replayMaxError = max(abs(replay.noTurbulence-reference.noTurbulence),[],'all');
        scatterMask = ~packets.ballistic;
        originalScatter = struct('pos',packets.pos(scatterMask,:), ...
            'dir',packets.dir(scatterMask,:), ...
            'powerW',packets.weight(scatterMask)*c.source.power/packets.photons, ...
            'pathLengthM',packets.pathLengthM(scatterMask), ...
            'sourcePhotons',packets.photons,'sourceSeed',reference.seed, ...
            'geometry',r.geometry,'referencePanel',[row,col]);
        save(fullfile(folder,'phase_scatter_original_r1_c1.mat'),'originalScatter','-v7.3');
        [X,Y] = ndgrid(r.x,r.y);
        ballisticField = bsf_phase_input(c,r.geometry,X,Y,zeros(size(X)));
        save(fullfile(folder,'phase_ballistic_r1_c1.mat'),'ballisticField','c');
        packetScatterPower = sum(originalScatter.powerW);
        referenceScatterPower = reference.powerWithoutTurbulence-reference.directPowerWithoutTurbulence;
    end
end
frequency = cell(2,1);
for j = 1:2
    row = panels(j,1);
    col = panels(j,2);
    o = loaded.s.panels{row,col}.out;
    [X,Y] = ndgrid(o.x(1:3:end),o.y(1:3:end));
    coarse = bsf_paper(c,o.geometry,X,Y,FrequencyNodes=2049);
    fine = bsf_paper(c,o.geometry,X,Y,FrequencyNodes=4097);
    frequency{j} = struct('row',row,'column',col,'pointGridStride',3, ...
        'relativeL2',norm(coarse.IBSF-fine.IBSF,'fro')/norm(fine.IBSF,'fro'), ...
        'scatterRelativeL2',norm(coarse.Is-fine.Is,'fro')/norm(fine.Is,'fro'), ...
        'minimumScatter',min(fine.Is,[],'all'));
end
tests = runtests(fullfile(root,'bsfPhysicsTest.m'));
testSummary = table(string({tests.Name}).',[tests.Passed].',[tests.Failed].', ...
    [tests.Incomplete].',[tests.Duration].', ...
    'VariableNames',{'Name','Passed','Failed','Incomplete','Seconds'});
writetable(testSummary,fullfile(folder,'tests.csv'));
fid = fopen(fullfile(folder,'tests.txt'),'w');
fprintf(fid,'%d tests; %d passed; %d failed; %d incomplete.\n',numel(tests), ...
    nnz([tests.Passed]),nnz([tests.Failed]),nnz([tests.Incomplete]));
fclose(fid);
validation = struct('independent',vertcat(records{:}),'paperFrequency',vertcat(frequency{:}), ...
    'originalReplayMaxError',replayMaxError,'originalScatterPacketPowerW',packetScatterPower, ...
    'originalScatterReferencePowerW',referenceScatterPower, ...
    'testPassed',all([tests.Passed]),'testCount',numel(tests), ...
    'pass',all(cellfun(@(v) v.compatible,records)) && replayMaxError==0 && all([tests.Passed]));
save(fullfile(folder,'validation.mat'),'validation');
fid = fopen(fullfile(folder,'validation.json'),'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(validation,PrettyPrint=true));
disp(validation);
end
