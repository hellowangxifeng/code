function intensity=rp_deposit_power(frequency,power,c,g,reference,qCenter,lambda,x,y)
%RP_DEPOSIT_POWER Conservatively map spectral-bin power through exact layered rays.
arguments
    frequency (:,2) double {mustBeReal,mustBeFinite}
    power (:,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
    c (1,1) struct
    g (1,1) struct
    reference (1,2) double {mustBeReal,mustBeFinite}
    qCenter (1,2) double {mustBeReal,mustBeFinite}
    lambda (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    x (1,:) double {mustBeReal,mustBeFinite}
    y (1,:) double {mustBeReal,mustBeFinite}
end
assert(size(frequency,1)==numel(power),'rp:SpectralPower','Frequency and power counts differ.');
assert(numel(x)>=2 && numel(y)>=2,'rp:ReceiverGrid','Receiver axes need at least two samples.');
dx=x(2)-x(1); dy=y(2)-y(1); nx=numel(x); ny=numel(y);
gridPower=zeros(ny*nx,1);
for start=1:100000:size(frequency,1)
    ids=start:min(start+99999,size(frequency,1));
    [points,propagating]=rp_carrier_landing(c,g,reference,lambda*frequency(ids,:)+qCenter);
    near=propagating & points(:,1)>=x(1)-dx & points(:,1)<=x(end)+dx & ...
        points(:,2)>=y(1)-dy & points(:,2)<=y(end)+dy;
    if ~any(near), continue; end
    points=points(near,:); weight=power(ids(near));
    gx=(points(:,1)-x(1))/dx+1; gy=(points(:,2)-y(1))/dy+1;
    ix=floor(gx); iy=floor(gy); fx=gx-ix; fy=gy-iy;
    for ox=0:1
        for oy=0:1
            xx=ix+ox; yy=iy+oy; wx=(1-fx)*(ox==0)+fx*(ox==1); wy=(1-fy)*(oy==0)+fy*(oy==1);
            valid=xx>=1 & xx<=nx & yy>=1 & yy<=ny;
            if ~any(valid), continue; end
            index=yy(valid)+(xx(valid)-1)*ny;
            gridPower=gridPower+accumarray(index,weight(valid).*wx(valid).*wy(valid),[ny*nx 1],@sum,0);
        end
    end
end
intensity=gridPower/(dx*dy);
end
