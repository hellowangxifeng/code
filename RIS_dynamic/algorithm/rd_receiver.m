function rx = rd_receiver(c,origin,center,Z,s,options)
%RD_RECEIVER One fixed physical receiver and its jitter / pointing samples.
%   rx = rd_receiver(c,origin,center,Z,samples) places a disc of radius
%   c.rx.apertureRadius at (center, Z) on the horizontal receiver plane.
%   Its nominal boresight r points from origin (above the frozen lower-face
%   centre by the glass thickness, as rp_aperture_channel) to the disc
%   centre; ALL links share this receiver. Samples follow
%   RIS_phase/rp_aperture_channel exactly:
%     offset_s = sigma_p * z_pos                      (horizontal, per axis)
%     e_s      = mu_theta + sigma_theta * z_ang,  psi_s = |e_s|
%     n_s      = cos(psi_s) r + sin(psi_s) (e_s1 e1 + e_s2 e2)/psi_s
%   with (e1,e2,r) the right-handed basis e1 = y x r / |.|, e2 = r x e1.
%   The receiver does not follow the waves; it is aligned once per condition.
if nargin<6, options = c.study; end
n = s.n;
rx = struct('center',reshape(center,1,2),'Z',Z,'origin',reshape(origin,1,3), ...
    'valid',all(isfinite(center)),'ray',nan(1,3),'e1',nan(1,3),'e2',nan(1,3), ...
    'offset',options.positionSigmaM*s.positionNormal,'normal',nan(n,3),'psi',zeros(n,1), ...
    'order',(1:n).','fovHalfAngleRad',options.fovHalfAngleRad, ...
    'apertureRadiusM',c.rx.apertureRadius,'fovOutageExact',NaN,'fovOutageSample',NaN);
errors = options.angleBiasRad+options.angleSigmaRad*s.angleNormal;
psi = vecnorm(errors,2,2);
rx.psi = psi;
rx.fovOutageSample = mean(psi>options.fovHalfAngleRad);
if options.angleSigmaRad>0
    rx.fovOutageExact = ncx2cdf((options.fovHalfAngleRad/options.angleSigmaRad)^2,2, ...
        sum(options.angleBiasRad.^2)/options.angleSigmaRad^2,'upper');
else
    rx.fovOutageExact = double(norm(options.angleBiasRad)>options.fovHalfAngleRad);
end
% Compact sample blocks (x bands, then y) let rd_capture prune packets per block.
[~,order] = sortrows([floor(rx.offset(:,1)/(2*c.rx.apertureRadius)),rx.offset(:,2)]);
rx.order = order;
if ~rx.valid, return; end
ray = [rx.center-rx.origin(1:2),Z-rx.origin(3)];
ray = ray/norm(ray);
e1 = cross([0 1 0],ray);
e1 = e1/norm(e1);
e2 = cross(ray,e1);
tangent = errors./max(psi,realmin);
rx.normal = cos(psi).*ray+sin(psi).*(tangent(:,1).*e1+tangent(:,2).*e2);
rx.ray = ray; rx.e1 = e1; rx.e2 = e2;
end
