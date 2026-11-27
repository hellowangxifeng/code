function summary = bsf_export_summary(folder)
%BSF_EXPORT_SUMMARY Write machine-readable metrics from completed panels.
arguments
    folder (1,1) string = fullfile(fileparts(mfilename('fullpath')),'results')
end
files = dir(fullfile(folder,'panel_r*_c*.mat'));
rows = cell(numel(files),1);
used = 0;
for j = 1:numel(files)
    loaded = load(fullfile(files(j).folder,files(j).name),'result');
    r = loaded.result;
    if ~isfield(r,'transportMetrics')
        continue
    end
    m = r.transportMetrics;
    p = r.paperMetrics;
    b = r.paperVsTransport;
    used = used+1;
    rows{used} = struct('row',r.row,'column',r.col,'incidenceDeg',r.geometry.thetaDeg, ...
        'originalPhotons',r.referencePhotons,'modelPhotons',r.photons, ...
        'precisionPassed',r.precisionPassed,'numericalSEpct',100*r.relativeSE, ...
        'originalNoisePct',100*m.referenceNoiseL2,'powerModelW',m.powerW, ...
        'powerOriginalW',m.referencePowerW,'powerModelSEW',r.powerSE,'powerErrorPct',100*m.powerRelative, ...
        'powerZ',r.powerZ,'rawL2pct',100*m.relativeL2,'rawL1pct',100*m.relativeL1, ...
        'peakErrorPct',100*m.peakRelative,'centroidErrorMM',1000*m.centroidErrorM, ...
        'amplitudeErrorPct',100*m.amplitudeL2,'chiPerBin',m.meanSquaredZ, ...
        'coverage95pct',100*m.coverage95,'eligibleBins',m.eligibleBins, ...
        'paperPowerErrorPct',100*p.powerRelative,'paperRawL2pct',100*p.relativeL2, ...
        'paperVsTransportL2pct',100*b.relativeL2,'paperVsTransportPowerPct',100*b.powerRelative, ...
        'ballisticQuadraturePct',100*r.ballisticQuadratureRelative, ...
        'truncatedPowerW',r.truncatedPower);
end
if used>0
    summary = struct2table(vertcat(rows{1:used}));
    summary = sortrows(summary,["row","column"]);
    writetable(summary,fullfile(folder,'comparison.csv'));
else
    summary = table();
end
end
