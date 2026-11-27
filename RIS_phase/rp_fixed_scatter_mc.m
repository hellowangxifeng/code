function samples=rp_fixed_scatter_mc(c,d,target)
%RP_FIXED_SCATTER_MC Stream more photons for the weak off-axis U-RIS baseline.
% Original 300k photons provide only ~27 packets in this remote 1 m window.
% Keep two independent half-runs to quantify receiver-statistic uncertainty.
rp_setup(); folder=fullfile(c.output,'fixed_receiver'); file=fullfile(folder,'scatter_samples.mat');
blocks=100; photonsPerBlock=1000000;
digest=java.security.MessageDigest.getInstance('SHA-256');
spec=struct('source',c.source,'water',c.water,'geometry',d.geometry, ...
    'mc',c.mc,'target',target,'blocks',blocks,'photonsPerBlock',photonsPerBlock);
digest.update(unicode2native([jsonencode(spec),fileread([mfilename('fullpath') '.m']), ...
    fileread(which('wris_mc')),fileread(which('rp_trace')),fileread(which('rp_illuminate'))],'UTF-8'));
signature=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
if isfile(file)
    old=load(file,'samples','signature');
    if strcmp(old.signature,signature), samples=old.samples; return; end
end
positions=cell(blocks,1); weights=positions; blockIds=positions; energy=zeros(blocks,1);
cc=c; cc.mc.photons=photonsPerBlock; cc.grid.bins=[32 32]; cc.grid.planeScale=1;
for k=1:blocks
    [~,packets]=wris_mc(cc,d.geometry,c.mc.seed+10000+k);
    use=~packets.ballistic;
    scatter=struct('bottomXY',(packets.pos(use,:)-d.geometry.pRis.')*d.geometry.RRis, ...
        'directionLocal',packets.dir(use,:)*d.geometry.RRis, ...
        'powerW',packets.weight(use)*c.source.power/photonsPerBlock, ...
        'ballistic',false(sum(use),1));
    scatter.bottomXY=scatter.bottomXY(:,1:2);
    illumination=rp_illuminate(c,d.geometry,scatter);
    r=rp_trace(c,d.geometry,illumination,d.bank{1}.phase);
    inside=r.valid & all(abs(r.landing(:,1:2)-target)<=.65,2);
    positions{k}=r.landing(inside,1:2); weights{k}=r.powerW(inside)/blocks;
    blockIds{k}=repmat(k,sum(inside),1); energy(k)=packets.energyResidual;
    if mod(k,5)==0
        fid=fopen(fullfile(folder,'scatter_progress.json'),'w');
        fprintf(fid,'%s',jsonencode(struct('completed',k,'total',blocks,'photons',k*photonsPerBlock))); fclose(fid);
        fprintf('Fixed U-RIS scatter: %d/%d blocks\n',k,blocks);
    end
end
samples=struct('xy',vertcat(positions{:}),'powerW',vertcat(weights{:}), ...
    'blockId',vertcat(blockIds{:}),'blocks',blocks,'photons',blocks*photonsPerBlock, ...
    'energyResidual',energy,'target',target,'signature',signature);
save(file,'samples','signature','-v7.3');
end
