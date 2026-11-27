function samples=rp_fixed_scatter_batch(c,d,target,batchId)
%RP_FIXED_SCATTER_BATCH Independent photon stream with its own restart cache.
c.mc.seed=c.mc.seed+100000*batchId;
c.output=fullfile(c.output,'fixed_receiver','mc_batches',sprintf('batch%02d',batchId));
folder=fullfile(c.output,'fixed_receiver');if ~isfolder(folder),mkdir(folder);end
samples=rp_fixed_scatter_mc(c,d,target);
end
