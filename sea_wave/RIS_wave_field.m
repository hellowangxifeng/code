function [elevation,gradient] = RIS_wave_field(motion,points,time)
%RIS_WAVE_FIELD Incident elevation and analytic x/y slopes at paired points.
% points is N-by-2 global [x,y]. time is elapsed seconds, scalar or N-vector.
% The same global phases drive motion, surfaces, profiles and diagnostics.
    arguments
        motion (1,1) struct
        points (:,2) double {mustBeReal,mustBeFinite,mustBeNonempty}
        time double {mustBeReal,mustBeFinite,mustBeNonnegative,mustBeNonempty}
    end
    if ~isscalar(time) && numel(time)~=size(points,1)
        error('RIS_wave_field:Time','time must be scalar or have one value per point.');
    end
    wave = motion.wave;
    absoluteTime = motion.parameters.startTime+time(:).';
    elevation = zeros(size(points,1),1);
    gradient = zeros(size(points));
    block = max(1,min(64,floor(2e6/size(points,1))));
    for first = 1:block:numel(wave.omega)
        ids = first:min(first+block-1,numel(wave.omega));
        phase = wave.omega(ids)*absoluteTime-wave.kx(ids)*points(:,1).'- ...
            wave.ky(ids)*points(:,2).';
        field = wave.coefficient(ids).*exp(1i*phase);
        elevation = elevation+real(sum(field,1)).';
        if nargout>1
            gradient = gradient+real([-1i*wave.kx(ids).';-1i*wave.ky(ids).']*field).';
        end
    end
end
