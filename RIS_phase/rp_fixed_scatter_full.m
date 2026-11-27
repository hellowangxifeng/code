function samples=rp_fixed_scatter_full(c,d,target)
%RP_FIXED_SCATTER_FULL Pool six independent 100-million-photon streams.
rp_setup();file=fullfile(c.output,'fixed_receiver','scatter_full.mat');
digest=java.security.MessageDigest.getInstance('SHA-256');
spec=struct('source',c.source,'water',c.water,'mc',c.mc,'geometry',d.geometry,'target',target);
digest.update(unicode2native([jsonencode(spec),fileread([mfilename('fullpath') '.m']), ...
    fileread(which('rp_fixed_scatter_batch')),fileread(which('rp_fixed_scatter_mc'))],'UTF-8'));
signature=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
if isfile(file)
    old=load(file,'samples','signature');
    if strcmp(old.signature,signature),samples=old.samples;return;end
end
parts=cell(1,6);parts{1}=rp_fixed_scatter_mc(c,d,target);
pool=gcp('nocreate');if isempty(pool),pool=parpool('Processes',3);end
jobs=parallel.FevalFuture.empty(0,5);
for k=1:5,jobs(k)=parfeval(pool,@rp_fixed_scatter_batch,1,c,d,target,k);end
for k=1:5,[id,value]=fetchNext(jobs);parts{id+1}=value;end
samples=parts{1};samples.xy=[];samples.powerW=[];samples.blockId=[];samples.energyResidual=[];
samples.sourceSignatures=cell(1,6);
for k=1:6
    samples.xy=[samples.xy;parts{k}.xy];samples.powerW=[samples.powerW;parts{k}.powerW/6];
    samples.blockId=[samples.blockId;parts{k}.blockId+(k-1)*100];
    samples.energyResidual=[samples.energyResidual;parts{k}.energyResidual];
    samples.sourceSignatures{k}=parts{k}.signature;
end
% Fix the density-estimation resolution from the pilot sample, independent of
% ABER/OP results. Increasing photons then improves variance at that resolution.
weights=parts{1}.powerW;normWeights=weights/sum(weights);xy=parts{1}.xy;
center=sum(xy.*normWeights,1);delta=xy-center;
covariance=(delta.*normWeights).'*delta;
effective=sum(weights)^2/sum(weights.^2);
samples.pilotBandwidthM=sqrt(trace(covariance)/2)*effective^(-1/6);
samples.blocks=600;samples.photons=600000000;samples.signature=signature;
save(file,'samples','signature','-v7.3');
end
