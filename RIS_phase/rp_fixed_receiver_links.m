function data=rp_fixed_receiver_links(d,c,policy)
%RP_FIXED_RECEIVER_LINKS Reproduce Fig2/3's common target and per-power search.
if nargin<3,policy='per_power';end
rp_setup();c.study.fovHalfAngleRad=.130;c.study.signature='fixed-receiver-v1';
signature=rp_fixed_receiver_signature(c,d,policy);
file=fullfile(c.output,'fixed_receiver','comparison.mat');
if isfile(file)
    old=load(file,'data');
    if isfield(old.data,'fixedReceiverSignature') && strcmp(old.data.fixedReceiverSignature,signature)
        data=old.data;return
    end
end
design=struct('distanceM',d.distanceM,'geometry',d.geometry,'illumination',d.illumination, ...
    'mc',d.mc,'bank',{{struct('phase',d.bank{1}.phase),struct('phase',d.bank{2}.phase)}});
pool=gcp('nocreate');if isempty(pool),pool=parpool('Processes',3);end
families={'steering'};widths=0;
for family={'flat_top','robust'}
    for width=[.01 .02 .03 .04:.02:.32]
        families{end+1}=family{1};widths(end+1)=width; %#ok<AGROW>
    end
end
jobs=parallel.FevalFuture.empty(0,numel(widths));
for k=1:numel(widths),jobs(k)=parfeval(pool,@rp_fixed_receiver_case,1,c,design,families{k},widths(k),.012);end
bank=cell(1,numel(widths));
for k=1:numel(widths),[id,b]=fetchNext(jobs);bank{id}=b;end
target=bank{1}.receiverCenterXY;
unmodulated=rp_fixed_receiver_unmodulated(c,design,target);
scatter=rp_fixed_scatter_full(c,design,target);
unmodulated=rp_fixed_refine_unmodulated(c,design,unmodulated,scatter);
while true
    data=rp_fixed_receiver_collect(c,design,bank,unmodulated,policy);
    selected=unique(data.choice(:));selected=selected(selected>0);pending=[];
    for k=selected.'
        if bank{k}.spot.r50M<.025 && ~isfield(bank{k},'refinement'),pending(end+1)=k;end %#ok<AGROW>
    end
    if isempty(pending),break;end
    for k=pending
        r=rp_fixed_receiver_refinement(c,design,bank{k});bank{k}=r.b;
    end
end
data.fixedReceiverSignature=signature;
data.sourceModelSignature=signature;
save(file,'data','-v7.3');save(fullfile(c.output,'fixed_receiver','bank.mat'),'bank','-v7.3');
end
