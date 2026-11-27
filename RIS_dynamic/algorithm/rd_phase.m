function p = rd_phase(c,g,a,base,family,width,target)
%RD_PHASE Realizable LC program of one RIS design, aimed at the common RX.
%   p = rd_phase(c,g,a,base,family,width,target), base = rp_phase(c,g,a).
%   Every tile m gets an affine phase (RIS_phase eq. 30 with H_m = 0)
%       phi_m(xi) = phi_0m + k0 b_m . (xi - r_m),
%   with b_m clamped and quantized per axis by rp_command (+-15.414 deg,
%   0.021 deg steps). The families differ only in WHERE each tile is sent:
%
%   'steering'  one global ramp b (rp_fixed_receiver_phase): the
%               power-weighted tile centroid is steered onto target, the
%               incident divergence is kept, phi_0m = k0 b.r_m (continuous).
%   'flat_top'  finite-width transport (rp_finite_phase): tile m is sent
%   'robust'    to target + eta_m, eta_m the inverse-CDF radius of its
%               whitened energy quantile u_m (flat: w sqrt(u); robust:
%               radial Weibull, p = c.study.robustExponent) at its own
%               azimuth. Small w converges all tiles onto the RX (focusing);
%               large w spreads them (jitter tolerance).
%               phi_0m = -Phi_in(r_m) - k0 L_m: equal optical path to EACH
%               tile's own landing point (RIS_phase eq. 50).
%
%   The only difference from RIS_phase is that the spot centre is the given
%   common receiver instead of each design's own reachable carrier.
k0 = 2*pi/c.source.wavelength;
p = base;
n = numel(p.unitIds);
ref = p.referenceXY;
weights = a.unitDirectPowerW(p.unitIds);
weights = weights/max(sum(weights),realmin);
target = reshape(target,1,2);
switch family
    case 'steering'
        centroid = sum(ref.*weights,1);
        incoming = rp_incident_field(c,g,centroid);
        outgoing = rp_paraxial(c,g,centroid,target);
        [b,command,saturated] = rp_command(c,outgoing.q-incoming.q);
        p.gradientOverK0 = repmat(b,n,1);
        p.commandDeg = repmat(command,n,1);
        saturated = repmat(saturated,n,1);
        p.phaseAtReferenceRad = mod(k0*(ref*b.'),2*pi);
        requested = repmat(target,n,1);
    case {'flat_top','robust'}
        requested = target+transportOffsets(c,ref,weights,family,width);
        wanted = rp_paraxial(c,g,ref,requested);
        [b,command,saturated] = rp_command(c,wanted.q-p.qInAtReference);
        p.gradientOverK0 = b;
        p.commandDeg = command;
    otherwise
        error('rd:Family','Unknown RIS design family "%s".',family);
end
p.phaseCurvature(:) = 0;
q = p.qInAtReference+p.gradientOverK0;
valid = sum(q.^2,2)<c.optics.nAir^2;
actual = requested;
for m = find(valid).'
    [actual(m,:),valid(m)] = rp_carrier_landing(c,g,ref(m,:),q(m,:));
end
actual(~valid,:) = requested(~valid,:);
p.targetsXY = actual;
p.outgoing = rp_paraxial(c,g,ref,actual);
p.outgoing.valid = p.outgoing.valid & valid;
if ~strcmp(family,'steering')
    input = rp_incident_field(c,g,ref);
    p.phaseAtReferenceRad = mod(-input.phaseRad-k0*p.outgoing.opl,2*pi);
end
p.phaseModel = 'linear';
p.family = family;
p.targetWidthM = width;
p.targetCenterXY = target;
p.requestedTargetsXY = requested;
p.saturated = saturated;
p.saturatedPowerFraction = sum(weights(saturated));
p.referenceDirectionWorld = [q,sqrt(max(0,c.optics.nAir^2-sum(q.^2,2)))]/c.optics.nAir*g.RRis.';
p.objective = sprintf('%s design aimed at the common receiver',family);
end

function offset = transportOffsets(c,ref,weights,family,width)
% RIS_phase rp_finite_phase: whitened energy quantiles -> target radii.
meanXY = sum(ref.*weights,1);
D = ref-meanXY;
covariance = (D.*weights).'*D;
[V,L] = eig((covariance+covariance.')/2);
W = V*diag(1./sqrt(max(diag(L),1e-12)))*V.';
white = D*W;
radius = vecnorm(white,2,2);
[~,order] = sort(radius);
cumulative = cumsum(weights(order));
u = zeros(numel(weights),1);
u(order) = cumulative-weights(order)/2;
u = max(1e-7,min(1-1e-7,u));
azimuth = atan2(white(:,2),white(:,1));
if strcmp(family,'flat_top')
    radial = width*sqrt(u);
else
    e = c.study.robustExponent;
    radial = width*(-log1p(-u)).^(1/e)/sqrt(2*gamma(1+2/e));
end
offset = radial.*[cos(azimuth),sin(azimuth)];
end
