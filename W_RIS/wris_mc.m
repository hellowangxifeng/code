function [out,packets] = wris_mc(c,g,seed,gain)
%WRIS_MC Generalized fig3-7 implicit-capture MC for a tilted finite detector.
% Reuses Table_abc through config and Hen_Green for the default HG g=0.924.
% Free path -log(U)/(a+b); at collision absorbed=w*a/c, remaining=w*b/c.
% Boundary is tested BEFORE collision so the exit segment has no fictitious
% absorption. No extra cos(theta) is applied: oblique ray landing density
% already supplies projected-area dilution. Each launch has power P0/N.
% The receiver is a first-passage scoring plane: off-board rays terminate in
% a separate ledger. Source-plane backscatter exits; outer path cap is logged.
if nargin<4, gain=1; end
old=rng; guard=onCleanup(@() rng(old));
rng(seed,'twister');
N=c.mc.photons; a=c.water.absorption; b=c.water.scattering; ext=a+b;
validateattributes(N,{'numeric'},{'scalar','integer','positive'});
validateattributes(c.mc.chunkSize,{'numeric'},{'scalar','integer','positive'});
validateattributes([a b],{'numeric'},{'real','finite','nonnegative'});
validateattributes(c.water.anisotropy,{'numeric'},{'scalar','>',-1,'<',1});
nx=c.grid.bins(1); ny=c.grid.bins(2);
xe=linspace(-c.ris.Lx/2,c.ris.Lx/2,nx+1); ye=linspace(-c.ris.Ly/2,c.ris.Ly/2,ny+1);
area=(xe(2)-xe(1))*(ye(2)-ye(1)); hist=zeros(nx,ny); hist2=hist;
% Same bin area on an extended scoring plane. Board acceptance and the
% energy ledger still use the ORIGINAL finite aperture, never this window.
planeScale=c.grid.planeScale;
validateattributes(planeScale,{'numeric'},{'scalar','integer','>=',1});
pxe=linspace(-planeScale*c.ris.Lx/2,planeScale*c.ris.Lx/2,planeScale*nx+1);
pye=linspace(-planeScale*c.ris.Ly/2,planeScale*c.ris.Ly/2,planeScale*ny+1);
planeHist=zeros(planeScale*nx,planeScale*ny);
ledger=struct('inside',0,'outside',0,'backEscape',0,'absorbed',0,'truncated',0,'rouletteDelta',0);
counts=struct('inside',0,'outside',0,'backEscape',0,'rouletteKilled',0,'truncated',0);
xyAll=cell(ceil(N/c.mc.chunkSize),1); wAll=xyAll;
xyDirect=xyAll; wDirect=xyAll; localDirectionSum=zeros(1,3);
xyCorePlane=xyAll; wCorePlane=xyAll;
positionPackets=xyAll; directionPackets=xyAll; weightPackets=xyAll;
ballisticPackets=xyAll; lengthPackets=xyAll;
for first=1:c.mc.chunkSize:N
    chunk=ceil(first/c.mc.chunkSize); n=min(c.mc.chunkSize,N-first+1);
    % Gaussian waist Wigner distribution: position and transverse slope are
    % INDEPENDENT. Correlating them launches an already diverging wavefront.
    uv=c.source.waist/2*randn(n,2);
    slope=c.source.divergence/2*randn(n,2);
    launch=[slope,ones(n,1)]; launch=launch./vecnorm(launch,2,2);
    d=launch*g.RTx.';
    p=[uv,zeros(n,1)]*g.RTx.'+g.pTx.';
    w=ones(n,1); path=zeros(n,1); scatter=zeros(n,1);
    active=true(n,1); hitXY=zeros(n,2); hitW=zeros(n,1); direct=false(n,1); accepted=false(n,1);
    hitLength=zeros(n,1);
    while any(active)
        id=find(active); pp=p(id,:); dd=d(id,:); ww=w(id);
        if ext>0, free=-log(max(rand(numel(id),1),realmin))/ext;
        else, free=inf(numel(id),1); end
        remaining=c.mc.maxPath-path(id); travel=min(free,remaining);
        [local,dist]=wris_intersect(pp,dd,g.pRis,g.RRis,[c.ris.Lx c.ris.Ly]);
        toward=dd*g.normal>1e-12;
        dist(~toward)=Inf;
        denom=dd*g.dTx;
        back=((g.pTx.'-pp)*g.dTx)./denom;
        back(denom>=-1e-12 | back<1e-10)=Inf;
        reach=dist<=travel & dist<=back & isfinite(dist);
        escape=back<=travel & back<dist & isfinite(back);
        at=local(reach,1:2); wr=ww(reach);
        inside=abs(at(:,1))<=c.ris.Lx/2 & abs(at(:,2))<=c.ris.Ly/2;
        rid=id(reach); captured=rid(inside);
        % Keep the first intersection BEFORE the finite-aperture test. Yaw
        % is measured from the complete core, including rays lost off-board.
        hitXY(rid,:)=at; hitW(rid)=wr;
        hitLength(rid)=path(rid)+dist(reach);
        direct(rid)=scatter(rid)==0; accepted(captured)=true;
        ledger.inside=ledger.inside+sum(wr(inside)); ledger.outside=ledger.outside+sum(wr(~inside));
        counts.inside=counts.inside+nnz(inside); counts.outside=counts.outside+nnz(~inside);
        if any(inside)
            localDirectionSum=localDirectionSum+sum((dd(reach,:)*g.RRis).* (wr.*inside),1);
        end
        ledger.backEscape=ledger.backEscape+sum(ww(escape));
        counts.backEscape=counts.backEscape+nnz(escape);
        active(id(reach|escape))=false;
        move=~reach & ~escape;
        if ~any(move), continue; end
        mid=id(move);
        p(mid,:)=pp(move,:)+dd(move,:).*travel(move); path(mid)=path(mid)+travel(move);
        cutoff=free(move)>=remaining(move); cut=mid(cutoff);
        ledger.truncated=ledger.truncated+sum(w(cut)); counts.truncated=counts.truncated+numel(cut);
        active(cut)=false; collide=mid(~cutoff);
        if isempty(collide), continue; end
        ledger.absorbed=ledger.absorbed+sum(w(collide))*a/ext;
        w(collide)=w(collide)*b/ext; scatter(collide)=scatter(collide)+1;
        low=collide(w(collide)<c.mc.rouletteThreshold & w(collide)>0);
        before=w(low); survive=rand(numel(low),1)<c.mc.rouletteSurvival;
        w(low)=before.*survive/c.mc.rouletteSurvival;
        ledger.rouletteDelta=ledger.rouletteDelta+sum(w(low)-before);
        counts.rouletteKilled=counts.rouletteKilled+nnz(~survive);
        active(collide(w(collide)==0))=false;
        cut=collide(scatter(collide)>=c.mc.maxScatter & active(collide));
        ledger.truncated=ledger.truncated+sum(w(cut)); counts.truncated=counts.truncated+numel(cut); active(cut)=false;
        sid=collide(active(collide)); ns=numel(sid);
        if ns==0, continue; end
        if c.water.anisotropy==0.924
            theta=Hen_Green(ns,2).'; % directly reuse original project's HG sampler
        else
            u=rand(ns,1); gg=c.water.anisotropy;
            if abs(gg)<1e-12, mu=2*u-1;
            else, mu=(1+gg^2-((1-gg^2)./(1-gg+2*gg*u)).^2)/(2*gg); end
            theta=acos(max(-1,min(1,mu)));
        end
        phi=2*pi*rand(ns,1); oldDir=d(sid,:);
        % Pole-safe orthonormal form of legacy direction-cosine update.
        helper=repmat([0 0 1],ns,1); helper(abs(oldDir(:,3))>0.9,:)=repmat([1 0 0],nnz(abs(oldDir(:,3))>0.9),1);
        e1=cross(helper,oldDir,2); e1=e1./vecnorm(e1,2,2); e2=cross(oldDir,e1,2);
        newDir=cos(theta).*oldDir+sin(theta).*(cos(phi).*e1+sin(phi).*e2);
        d(sid,:)=newDir./vecnorm(newDir,2,2);
    end
    valid=accepted & hitW>0; xyAll{chunk}=hitXY(valid,:); wAll{chunk}=hitW(valid);
    % Bin each packet once per chunk, after transport; no dense histogram
    % allocation at every scattering step is needed.
    if any(valid)
        ix=discretize(hitXY(valid,1),xe); iy=discretize(hitXY(valid,2),ye);
        hist=hist+accumarray([ix iy],hitW(valid),[nx ny]);
        hist2=hist2+accumarray([ix iy],hitW(valid).^2,[nx ny]);
    end
    landed=hitW>0;
    pix=discretize(hitXY(landed,1),pxe); piy=discretize(hitXY(landed,2),pye);
    shown=isfinite(pix) & isfinite(piy); landedW=hitW(landed);
    if any(shown)
        planeHist=planeHist+accumarray([pix(shown) piy(shown)],landedW(shown),size(planeHist));
    end
    core=accepted & direct;
    xyDirect{chunk}=hitXY(core,:); wDirect{chunk}=hitW(core);
    xyCorePlane{chunk}=hitXY(direct,:); wCorePlane{chunk}=hitW(direct);
    if nargout>1
        positionPackets{chunk}=[hitXY(valid,:),zeros(nnz(valid),1)]*g.RRis.'+g.pRis.';
        directionPackets{chunk}=d(valid,:); weightPackets{chunk}=hitW(valid);
        ballisticPackets{chunk}=direct(valid); lengthPackets{chunk}=hitLength(valid);
    end
end
xy=vertcat(xyAll{:}); weights=vertcat(wAll{:});
scale=c.source.power/N;
out.irradiance=hist*scale/area*gain; out.noTurbulence=hist*scale/area;
out.standardError=sqrt(max(0,hist2-hist.^2/N)/(N*(N-1)))*c.source.power/area*gain;
out.x=(xe(1:end-1)+xe(2:end))/2; out.y=(ye(1:end-1)+ye(2:end))/2;
out.cellArea=area; out.capturedPower=ledger.inside*scale*gain;
out.plane=struct('x',(pxe(1:end-1)+pxe(2:end))/2, ...
    'y',(pye(1:end-1)+pye(2:end))/2,'noTurbulence',planeHist*scale/area, ...
    'powerWithoutTurbulence',(ledger.inside+ledger.outside)*scale);
out.powerWithoutTurbulence=ledger.inside*scale;
out.powerSE=c.source.power*sqrt(max(0,sum(weights.^2)-sum(weights)^2/N)/(N*(N-1)));
out.peak=max(out.irradiance,[],'all'); out.gain=gain;
out.metrics=wris_metrics(xy,weights);
out.directMetrics=wris_metrics(vertcat(xyDirect{:}),vertcat(wDirect{:}));
out.directPlaneMetrics=wris_metrics(vertcat(xyCorePlane{:}),vertcat(wCorePlane{:}));
out.directPowerWithoutTurbulence=sum(vertcat(wDirect{:}))*scale;
out.directPlanePowerWithoutTurbulence=sum(vertcat(wCorePlane{:}))*scale;
out.coreClippedFraction=1-out.directPowerWithoutTurbulence/max(out.directPlanePowerWithoutTurbulence,realmin);
out.meanIncidentDirectionLocal=localDirectionSum/max(ledger.inside,realmin);
out.ledger=ledger; out.counts=counts; out.seed=seed; out.photons=N;
out.energyResidual=(ledger.inside+ledger.outside+ledger.backEscape+ledger.absorbed+ledger.truncated-N-ledger.rouletteDelta)/N;
out.geometry=g;
if nargout>1
    % Packet path length is transport bookkeeping, not a coherent scattering phase.
    packets=struct('pos',vertcat(positionPackets{:}),'dir',vertcat(directionPackets{:}), ...
        'weight',vertcat(weightPackets{:}),'ballistic',vertcat(ballisticPackets{:}), ...
        'pathLengthM',vertcat(lengthPackets{:}),'photons',N, ...
        'energyResidual',out.energyResidual,'ledger',ledger);
end
end
