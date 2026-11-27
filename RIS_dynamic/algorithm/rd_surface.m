function sea = rd_surface(motion,time,settings,origin)
%RD_SURFACE Frozen Fourier sea surface: cubic grid near the beam, exact sum elsewhere.
wave = motion.wave;
a = wave.coefficient(:).*exp(1i*wave.omega(:)*(motion.parameters.startTime+time));
k = [wave.kx(:),wave.ky(:)];
active = abs(a)>0;
a = a(active);
a = a(:);
k = k(active,:);
sea.amplitudeBound = sum(abs(a));
sea.slopeBound = (2*sum(abs(a).*abs(k),1)); % conservative for cubic interpolation
sea.gridBounds = [];
sea.evaluate = @(xy) evaluate(xy,a,k,[],[]);
if isempty(a) || settings.gridHalfWidth == 0
    return
end

n = ceil(settings.gridHalfWidth/settings.gridStep);
x = origin(1) + (-n:n)*settings.gridStep;
y = origin(2) + (-n:n)*settings.gridStep;
% Separable evaluation of the exact Fourier series at every grid point.
field = zeros(numel(x),numel(y));
for first = 1:256:numel(a)
    ids = first:min(first+255,numel(a));
    ex = exp(-1i*x(:)*k(ids,1).');
    ey = exp(-1i*k(ids,2)*y);
    field = field + real((ex.*a(ids).')*ey);
end
interpolant = griddedInterpolant({x,y},field,'cubic','none');
sea.gridBounds = [x(2),x(end-1);y(2),y(end-1)];
sea.evaluate = @(xy) evaluate(xy,a,k,interpolant,sea.gridBounds);
% Cubic convolution can overshoot node extrema. These factors bound the
% interior cubic-convolution weights and derivative of adjacent differences.
sea.gridAmplitudeBound = 2*max(abs(field),[],'all');
sea.gridSlopeBound = 3*[max(abs(diff(field,1,1)),[],'all'), ...
    max(abs(diff(field,1,2)),[],'all')]/settings.gridStep;
sea.amplitudeBound = sea.amplitudeBound + settings.gridStep*norm(sea.slopeBound);
end

function [eta,gradient] = evaluate(xy,a,k,interpolant,interior)
eta = zeros(size(xy,1),1);
gradient = zeros(size(xy));
if isempty(a), return; end
inside = false(size(eta));
if ~isempty(interpolant)
    inside = all(xy>=interior(:,1).' & xy<=interior(:,2).',2);
    x = xy(inside,1);
    y = xy(inside,2);
    eta(inside) = interpolant(x,y);
    if nargout > 1
        epsM = 1e-4;
        gradient(inside,:) = [(interpolant(x+epsM,y)-interpolant(x-epsM,y)), ...
            (interpolant(x,y+epsM)-interpolant(x,y-epsM))]/(2*epsM);
    end
end
outside = find(~inside);
for first = 1:128:numel(outside)
    ids = outside(first:min(first+127,numel(outside)));
    field = exp(-1i*xy(ids,:)*k.').*a.';
    eta(ids) = real(sum(field,2));
    if nargout > 1, gradient(ids,:) = imag(field*k); end
end
end
