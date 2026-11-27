function [intensity,errorBound] = rp_pixel_intensity(fields,slopeX,slopeY,pixelSize,maxError)
%RP_PIXEL_INTENSITY Pixel-integrated coherent superposition of local plane carriers.
arguments
    fields (:,:) {mustBeNumeric}
    slopeX (:,:) double {mustBeReal,mustBeFinite}
    slopeY (:,:) double {mustBeReal,mustBeFinite}
    pixelSize (1,2) double {mustBeReal,mustBePositive,mustBeFinite}
    maxError (1,1) double {mustBeReal,mustBeNonnegative,mustBeFinite} = 0
end
assert(isequal(size(fields),size(slopeX),size(slopeY)),'rp:CoherenceSize','Carrier arrays must match fields.');
% Empty receiver-window columns contribute exactly zero. Removing them does
% not approximate coherence and avoids quadratic work on out-of-window tiles.
active=any(fields~=0,1);
fields=fields(:,active); slopeX=slopeX(:,active); slopeY=slopeY(:,active);
intensity=sum(abs(fields).^2,2);
first=zeros(1,size(fields,2)); last=first;
for m=1:size(fields,2)
    support=find(fields(:,m)~=0);
    first(m)=support(1); last(m)=support(end);
end
% Optional, rigorously bounded removal of negligible cross terms. A common
% carrier is subtracted only to tighten slope intervals; differences are unchanged.
count=size(fields,2); skip=false(count); errorBound=0;
if maxError>0 && count>1
    [~,reference]=max(sum(abs(fields).^2,1));
    commonX=slopeX(:,reference); commonY=slopeY(:,reference);
    lowX=zeros(count,1); highX=lowX; lowY=lowX; highY=lowX;
    peak=max(abs(fields),[],1).';
    for m=1:count
        present=fields(:,m)~=0;
        x=slopeX(present,m)-commonX(present); y=slopeY(present,m)-commonY(present);
        lowX(m)=min(x); highX(m)=max(x); lowY(m)=min(y); highY(m)=max(y);
    end
    gapX=max(0,max(lowX-highX.',lowX.'-highX));
    gapY=max(0,max(lowY-highY.',lowY.'-highY));
    factor=min(1,2./max(gapX*pixelSize(1),realmin)).* ...
        min(1,2./max(gapY*pixelSize(2),realmin));
    bound=2*(peak*peak.').*factor;
    [row,column]=find(triu(true(count),1));
    linear=sub2ind([count count],row,column);
    [value,order]=sort(bound(linear));
    selected=order(cumsum(value)<=maxError);
    skip(linear(selected))=true;
    errorBound=sum(value(1:numel(selected)));
end
for m=1:size(fields,2)-1
    for n=m+1:size(fields,2)
        if skip(m,n) || last(m)<first(n) || last(n)<first(m), continue; end
        overlap=fields(:,m)~=0 & fields(:,n)~=0;
        if ~any(overlap), continue; end
        dx=(slopeX(overlap,m)-slopeX(overlap,n))*pixelSize(1)/(2*pi);
        dy=(slopeY(overlap,m)-slopeY(overlap,n))*pixelSize(2)/(2*pi);
        factor=normalisedSinc(dx).*normalisedSinc(dy);
        intensity(overlap)=intensity(overlap)+2*real(fields(overlap,m).*conj(fields(overlap,n))).*factor;
    end
end
intensity=max(real(intensity),0);
end

function value=normalisedSinc(x)
value=ones(size(x)); nz=abs(x)>1e-12;
value(nz)=sin(pi*x(nz))./(pi*x(nz));
end
