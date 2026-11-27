function transmission = rp_fresnel(nFrom,nTo,q)
%RP_FRESNEL Unpolarised Fresnel power transmission, zero for TIR/invalid modes.
arguments
    nFrom (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    nTo (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
    q (:,2) double {mustBeReal,mustBeFinite}
end
q2=sum(q.^2,2); k1=sqrt(max(0,nFrom^2-q2)); k2=sqrt(max(0,nTo^2-q2));
rs=((k1-k2)./max(k1+k2,realmin)).^2;
rp=((nTo^2*k1-nFrom^2*k2)./max(nTo^2*k1+nFrom^2*k2,realmin)).^2;
transmission=max(0,min(1,1-(rs+rp)/2));
transmission(q2>=min(nFrom^2,nTo^2))=0;
end
