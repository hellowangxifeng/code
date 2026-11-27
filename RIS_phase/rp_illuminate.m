function a = rp_illuminate(c,g,mc)
%RP_ILLUMINATE Apply physical LC apertures, retaining Gaussian tails and scattering.
arguments
    c (1,1) struct
    g (1,1) struct
    mc (1,1) struct
end
beam=bottomEllipse(c,g);
q=c.optics.nWater*mc.directionLocal(:,1:2);
kz=sqrt(c.optics.nGlass^2-sum(q.^2,2));
lc=mc.bottomXY+g.hIn*q./kz;
count=c.array.count; pitch=c.array.pitchM;
ix=round(lc(:,1)/pitch(1)+(count(1)+1)/2);
iy=round(lc(:,2)/pitch(2)+(count(2)+1)/2);
center=[(ix-(count(1)+1)/2)*pitch(1),(iy-(count(2)+1)/2)*pitch(2)];
inside=ix>=1 & ix<=count(1) & iy>=1 & iy<=count(2) & ...
    all(abs(lc-center)<=c.array.apertureM/2,2) & ...
    all(abs(lc)<=[c.ris.Lx c.ris.Ly]/2,2);
id=ix(inside)+(iy(inside)-1)*count(1);
[ux,uy]=meshgrid(((1:count(1))-(count(1)+1)/2)*pitch(1), ...
    ((1:count(2))-(count(2)+1)/2)*pitch(2));
centers=[reshape(ux.',[],1),reshape(uy.',[],1)];
direct=mc.ballistic(inside);
power=mc.powerW(inside).*rp_fresnel(c.optics.nWater,c.optics.nGlass,q(inside,:));
unitPower=accumarray(id,power,[prod(count) 1],@sum,0);
directPower=accumarray(id(direct),power(direct),[prod(count) 1],@sum,0);
ids=(1:prod(count)).'; reference=centers;
for dimension=1:2
    moment=accumarray(id,power.*lc(inside,dimension),[prod(count) 1],@sum,0);
    use=unitPower>0;
    reference(use,dimension)=moment(use)./unitPower(use);
end
lookup=zeros(prod(count),1); lookup(ids)=1:numel(ids);
% PDF Eq.(8) supplies deterministic ballistic power; MC remains an independent
% estimator for scattered paths. No fitted optical-power normalization.
unitDirectMC=directPower;
[nodes,weights]=rp_gauss_legendre(c.input.apertureOrder);
[qx,qy]=meshgrid(nodes*c.array.apertureM(1)/2,nodes*c.array.apertureM(2)/2);
quadWeights=weights*weights.'/4;
quad=repelem(centers(ids,:),numel(qx),1)+repmat([qx(:),qy(:)],numel(ids),1);
modelWeight=reshape(rp_bsf_ballistic(c,g,quad,g.hIn),numel(qx),numel(ids));
modelWeight=modelWeight.*quadWeights(:);
modelPower=sum(modelWeight,1).'*prod(c.array.apertureM);
scale=1;
directPower(:)=0; directPower(ids)=modelPower;
unitPower=unitPower-unitDirectMC+directPower;
for dimension=1:2
    coordinate=reshape(quad(:,dimension),numel(qx),numel(ids));
    centroid=sum(modelWeight.*coordinate,1).'./max(sum(modelWeight,1).',realmin);
    use=modelPower>max(modelPower)*1e-10;
    reference(use,dimension)=centroid(use);
end
% Match ballistic ray quadrature to each deterministic tile integral. Scatter
% packet weights retain their original P0/N and are never divided by N again.
packetScale=ones(size(power));
hasCore=direct & unitDirectMC(id)>0;
packetScale(hasCore)=directPower(id(hasCore))./unitDirectMC(id(hasCore));
power=power.*packetScale;
corners=[-1 -1;1 -1;1 1;-1 1;-1 -1].*c.array.apertureM/2;
points=repelem(centers,5,1)+repmat(corners,size(centers,1),1);
path=rp_layer_path(points-g.sourceLocal(1:2).',[-g.sourceLocal(3),g.hIn], ...
    [c.optics.nWater,c.optics.nGlass],c.numerics.bisectionIterations);
lower=g.sourceLocal(1:2).'-g.sourceLocal(3)*path.q./path.kz(:,1);
a=struct('ellipse',beam,'bottomXY',mc.bottomXY(inside,:),'lcXY',lc(inside,:), ...
    'qIn',q(inside,:),'powerW',power,'ballistic',direct,'unitId',id, ...
    'sampleUnitIndex',lookup(id),'unitCenters',centers,'activeIds',ids, ...
    'referenceXY',reference,'unitPowerW',unitPower,'unitDirectPowerW',directPower, ...
    'activePowerW',sum(unitPower),'activeDirectPowerW',sum(directPower), ...
    'apertureBottomX',reshape(lower(:,1),5,[]),'apertureBottomY',reshape(lower(:,2),5,[]));
a.unitDirectPowerMCW=unitDirectMC;
a.ballisticModelScale=scale;
a.incidentFieldModel='BSF Eq.(8) deterministic amplitude, assumed Gaussian phase; MC scatter';
end

function beam=bottomEllipse(c,g)
slope=tan(c.source.divergence);
cutoff=slope*sqrt(-log(c.source.edgeRelativeIntensity)/2);
cosSquared=1/(1+cutoff^2); t=g.directionLocal(1:2); depth=-g.sourceLocal(3);
matrix=cosSquared*eye(2)-t*t.';
assert(min(eig(matrix))>0,'rp:BeamEllipse','The lower-face cone intersection must be an ellipse.');
shift=matrix\(depth*g.directionLocal(3)*t);
radiusSquared=shift.'*matrix*shift-depth^2*(cosSquared-g.directionLocal(3)^2);
shape=radiusSquared*(matrix\eye(2));
[vectors,values]=eig(shape); [axesM,order]=sort(sqrt(diag(values)),'descend');
vectors=vectors(:,order); angle=linspace(0,2*pi,241);
center=g.sourceLocal(1:2).'+shift.';
boundary=(vectors*(axesM.*[cos(angle);sin(angle)])).'+center;
beam=struct('center',center,'shape',shape,'semiAxesM',axesM.', ...
    'boundaryXY',boundary,'axesVectors',vectors,'waistAtSourceM',c.source.waist, ...
    'beamRadiusAtRISM',norm(g.sourceLocal)*tan(c.source.divergence));
end
