function bsf_run_all(options)
%BSF_RUN_ALL Reproducible adaptive precision study, all 18 ORIGINAL snapshots.
% Each batch is saved; a repeated call resumes its exact independent seeds.
arguments
    options.TargetRelativeSE (1,1) {mustBeNumeric,mustBePositive,mustBeReal} = 0.001
    options.BatchPhotons (1,1) {mustBeNumeric,mustBeInteger,mustBePositive,mustBeReal} = 8e6
    options.MaxPhotons (1,1) {mustBeNumeric,mustBeInteger,mustBePositive,mustBeReal} = 1024e6
end
root = bsf_setup();
output = fullfile(root,'results');
if ~isfolder(output)
    mkdir(output);
end
source = fullfile(fileparts(root),'results','wris_results.mat');
loaded = load(source,'s');
s = loaded.s;
c = s.config;
versionTag = 'bsf-rte-v1-forced-b-continuous-a';
sourceInfo = dir(source);
save(fullfile(output,'input_snapshot.mat'),'c','source','sourceInfo','versionTag','options');
diary(fullfile(output,'run.log'));
cleanup = onCleanup(@() diary('off'));
fprintf('Start %s | target %.4f%% | original source %s\n',string(datetime('now')), ...
    100*options.TargetRelativeSE,source);
for row = 1:size(s.panels,1)
    for col = 1:size(s.panels,2)
        original = s.panels{row,col}.out;
        g = original.geometry;
        targetFile = fullfile(output,sprintf('panel_r%d_c%d.mat',row,col));
        if isfile(targetFile)
            saved = load(targetFile,'result');
            result = saved.result;
            assert(strcmp(result.versionTag,versionTag) && isequaln(result.config,c), ...
                'bsf:StaleCache','Saved model or configuration changed; use a separate results folder.');
        else
            result = initialise(c,g,original,versionTag,row,col);
        end
        fprintf('r%d c%d | theta %.4f deg | original power %.9g W\n',row,col,g.thetaDeg, ...
            original.powerWithoutTurbulence);
        while result.photons<options.MaxPhotons && result.relativeSE>options.TargetRelativeSE
            batch = result.batches+1;
            count = min(options.BatchPhotons,options.MaxPhotons-result.photons);
            seed = 9000000+100000*row+10000*col+batch;
            trial = bsf_scatter_transport(c,g,Photons=count,Seed=seed);
            result = combine(result,trial);
            result.targetRelativeSE = options.TargetRelativeSE;
            result.precisionPassed = result.relativeSE<=options.TargetRelativeSE;
            result.updatedAt = string(datetime('now'));
            save(targetFile,'result','-v7.3');
            writeProgress(output,result);
            fprintf('  N=%g | numerical map SE %.5f%% | power %.9g W | truncation %.3g W\n', ...
                result.photons,100*result.relativeSE,sum(result.IBSF,'all')*original.cellArea, ...
                result.truncatedPower);
        end
        result = finishMetrics(result,original);
        save(targetFile,'result','-v7.3');
        bsf_export_summary(output);
    end
end
progress = struct('status','complete','time',string(datetime('now')));
writeJson(fullfile(output,'progress.json'),progress);
bsf_plot_results(output);
fprintf('Completed %s\n',string(datetime('now')));
end

function result = initialise(c,g,o,versionTag,row,col)
Ib = bsf_cell_average(@(x,y) bsf_ballistic(c,g,x,y),o.x,o.y,4);
IbFine = bsf_cell_average(@(x,y) bsf_ballistic(c,g,x,y),o.x,o.y,6);
paper = bsf_paper_map(c,g,o.x,o.y);
result = struct('versionTag',versionTag,'config',c,'geometry',g,'row',row,'col',col, ...
    'x',o.x,'y',o.y,'cellArea',o.cellArea,'Ib',IbFine,'paper',paper, ...
    'ballisticQuadratureRelative',norm(Ib-IbFine,'fro')/max(norm(IbFine,'fro'),realmin), ...
    'photons',0,'batches',0,'weightedMean',zeros(size(Ib)), ...
    'm2',zeros(size(Ib)),'powerMean',0,'powerM2',0, ...
    'truncatedPower',0,'relativeSE',inf,'Is',zeros(size(Ib)), ...
    'IBSF',IbFine,'standardError',zeros(size(Ib)),'powerSE',0, ...
    'targetRelativeSE',NaN,'precisionPassed',false,'updatedAt',"", ...
    'referencePhotons',o.photons,'referenceSeed',o.seed,'referenceGain',o.gain);
end

function result = combine(result,trial)
% Exact pooled first and second moments for independent, not identically
% forced geometries: all samples HERE share the same frozen geometry.
n = result.photons;
m = trial.photons;
total = n+m;
delta = trial.Is-result.weightedMean;
result.m2 = result.m2+trial.standardError.^2*m*(m-1)+delta.^2*n*m/total;
result.weightedMean = result.weightedMean+delta*m/total;
deltaPower = trial.power-result.powerMean;
result.powerM2 = result.powerM2+trial.powerSE^2*m*(m-1)+deltaPower^2*n*m/total;
result.powerMean = result.powerMean+deltaPower*m/total;
result.truncatedPower = (result.truncatedPower*n+trial.truncatedPower*m)/total;
result.photons = total;
result.batches = result.batches+1;
result.Is = result.weightedMean;
result.IBSF = result.Ib+result.Is;
result.standardError = sqrt(result.m2/(total*(total-1)));
result.powerSE = sqrt(result.powerM2/(total*(total-1)));
result.relativeSE = norm(result.standardError,'fro')/norm(result.IBSF,'fro');
end

function result = finishMetrics(result,o)
referenceSE = o.standardError/o.gain;
result.transportMetrics = bsf_compare(result.IBSF,o.noTurbulence,o.x,o.y,o.cellArea, ...
    result.standardError,referenceSE);
result.paperMetrics = bsf_compare(result.paper.IBSF,o.noTurbulence,o.x,o.y,o.cellArea, ...
    zeros(size(result.Ib)),referenceSE);
result.paperVsTransport = bsf_compare(result.paper.IBSF,result.IBSF,o.x,o.y,o.cellArea, ...
    zeros(size(result.Ib)),result.standardError);
result.powerZ = (sum(result.IBSF,'all')*o.cellArea-o.powerWithoutTurbulence)/ ...
    hypot(result.powerSE,o.powerSE);
result.coreMasks = [0.1,0.25,0.5];
[X,Y] = ndgrid(o.x,o.y);
distance = hypot(X-result.geometry.axisLocal(1),Y-result.geometry.axisLocal(2));
result.encircledPowerW = zeros(3,3);
for j = 1:3
    mask = distance<=result.coreMasks(j);
    result.encircledPowerW(j,:) = [sum(result.IBSF(mask)),sum(o.noTurbulence(mask)), ...
        sum(result.paper.IBSF(mask))]*o.cellArea;
end
end

function writeProgress(folder,result)
progress = struct('status','running','row',result.row,'column',result.col, ...
    'photons',result.photons,'relativeNumericalSE',result.relativeSE, ...
    'target',result.targetRelativeSE,'time',result.updatedAt);
writeJson(fullfile(folder,'progress.json'),progress);
end

function writeJson(file,value)
fid = fopen(file,'w');
assert(fid>=0,'bsf:Output','Cannot open progress file.');
guard = onCleanup(@() fclose(fid));
fprintf(fid,'%s',jsonencode(value,PrettyPrint=true));
end
