function mc = rp_mc(c,g)
%RP_MC Reuse underwater absorption/scattering MC and score the lower-face map.
arguments
    c (1,1) struct
    g (1,1) struct
end
key=struct('source',c.source,'water',c.water,'mc',c.mc,'pTx',g.pTx,'RTx',g.RTx, ...
    'pBottom',g.pRis,'RRis',g.RRis,'size',[c.ris.Lx c.ris.Ly]);
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native(jsonencode(key),'UTF-8'));
digest.update(unicode2native(fileread(mfilename('fullpath')+".m"),'UTF-8'));
digest.update(unicode2native(fileread(which('wris_mc')),'UTF-8'));
digest.update(unicode2native(fileread(which('rp_incident_field')),'UTF-8'));
digest.update(unicode2native(fileread(which('rp_bsf_ballistic')),'UTF-8'));
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
cacheFile=fullfile(c.cache.folder,['mc_',hash,'.mat']);
if c.cache.enabled && isfile(cacheFile)
    stored=load(cacheFile,'mc'); mc=stored.mc; return
end
transport=c; transport.grid.bins=c.mc.mapBins; transport.grid.planeScale=1;
[~,packets]=wris_mc(transport,g,c.mc.seed);
local=(packets.pos-g.pRis.')*g.RRis;
direction=packets.dir*g.RRis;
xe=linspace(-c.ris.Lx/2,c.ris.Lx/2,c.mc.mapBins(1)+1);
ye=linspace(-c.ris.Ly/2,c.ris.Ly/2,c.mc.mapBins(2)+1);
ix=discretize(local(:,1),xe); iy=discretize(local(:,2),ye);
power=packets.weight*c.source.power/packets.photons;
area=(xe(2)-xe(1))*(ye(2)-ye(1));
map=accumarray([iy,ix],power,fliplr(c.mc.mapBins),@sum,0)/area;
core=packets.ballistic;
direct=accumarray([iy(core),ix(core)],power(core),fliplr(c.mc.mapBins),@sum,0)/area;
corePower=power(core); coreXY=local(core,1:2);
meanXY=sum(coreXY.*corePower,1)/max(sum(corePower),realmin);
deviation=coreXY-meanXY;
covariance=(deviation.*corePower).'*deviation/max(sum(corePower),realmin);
mc=struct('bottomXY',local(:,1:2),'directionLocal',direction,'powerW',power, ...
    'ballistic',core,'x',(xe(1:end-1)+xe(2:end))/2,'y',(ye(1:end-1)+ye(2:end))/2, ...
    'irradianceWm2',map,'ballisticIrradianceWm2',direct,'cellAreaM2',area, ...
    'incidentPowerW',sum(power),'photons',packets.photons, ...
    'energyResidual',packets.energyResidual,'ledger',packets.ledger, ...
    'ballisticMeanXY',meanXY,'ballisticCovarianceM2',covariance);
mc.pathLengthM=packets.pathLengthM;
[xx,yy]=meshgrid(mc.x,mc.y);
input=rp_incident_field(c,g,[xx(:),yy(:)],0);
normalisation=1;
mc.ballisticModelScale=normalisation;
mc.ballisticField=reshape(input.field*sqrt(normalisation),size(xx));
mc.ballisticPhaseRad=reshape(input.phaseRad,size(xx));
mc.ballisticModelIrradianceWm2=abs(mc.ballisticField).^2;
mc.incidentModelPowerW=sum(input.irradianceWm2)*area+sum(power(~core));
mc.fieldModel='BSF Eq.(8) amplitude; explicit Gaussian phase; incoherent HG scatter';
if c.cache.enabled
    if ~isfolder(c.cache.folder), mkdir(c.cache.folder); end
    save(cacheFile,'mc','-v7.3');
end
end
