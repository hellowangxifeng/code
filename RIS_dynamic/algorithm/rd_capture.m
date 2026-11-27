function gain = rd_capture(c,pk,rx)
%RD_CAPTURE Exact packet capture by a tilted finite disc with a per-ray FOV.
%   gain = rd_capture(c,pk,rx) returns, for every receiver sample s, the
%   captured power / P0. Packet j lands at X_j on the plane z = rx.Z with
%   unit direction d_j and power w_j. Sample s has disc centre P_s = C + delta_s
%   (on the same plane) and unit normal n_s. The packet line meets the disc
%   plane at Y = X_j + t d_j,  t = n_s.(P_s - X_j) / (n_s.d_j), and
%     captured_s = sum_j w_j 1{|Y - P_s| <= a} 1{d_j.n_s >= cos(Theta)}.
%   This is the packet form of RIS_phase's tilted-disc integral, but the FOV
%   test is applied to EVERY ray: after a rough sea the arrival directions
%   spread, so a chief-ray FOV test would be wrong for the direct link.
n = size(rx.offset,1);
gain = zeros(n,1);
if ~rx.valid || isempty(pk.power), return; end
a = rx.apertureRadiusM;
cosFov = cos(rx.fovHalfAngleRad);
x = pk.xy-rx.center;
d = pk.dir;
w = pk.power;
% A disc point Y lies within a of P and within a*sin(tilt) <= a of the plane,
% so its packet lands within a(1 + tan(theta_j)) of P horizontally.
reach = a*(1+min(vecnorm(d(:,1:2),2,2)./max(d(:,3),realmin),20));
cone = cos(min(pi/2,rx.fovHalfAngleRad+max(rx.psi)));
keep = d(:,3)>0 & d*rx.ray(:)>=cone-1e-12 & ...
    vecnorm(x,2,2)<=max(vecnorm(rx.offset,2,2))+reach;
if ~any(keep), return; end
x = x(keep,:); d = d(keep,:); w = w(keep); reach = reach(keep);
block = max(16,min(n,floor(c.capture.maxPairs/numel(w))));
for first = 1:block:n
    ids = rx.order(first:min(first+block-1,n));
    P = rx.offset(ids,:);
    N = rx.normal(ids,:);
    near = all(x>=min(P,[],1)-reach & x<=max(P,[],1)+reach,2);
    if ~any(near), continue; end
    X = x(near,:); D = d(near,:); W = w(near);
    dn = N*D.';                                            % n_s . d_j
    t = (sum(N(:,1:2).*P,2)-N(:,1:2)*X.')./dn;             % to the disc plane
    r2 = sum(X.^2,2).'+sum(P.^2,2)-2*(P*X.') + ...
        2*t.*(sum(X.*D(:,1:2),2).'-P*D(:,1:2).')+t.^2;     % |Y - P_s|^2
    hit = dn>=cosFov & r2<=a^2;
    gain(ids) = double(hit)*W;
end
gain = gain/c.source.power;
end
