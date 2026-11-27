function [X,Y,Z] = RIS_sea_surface(motion,x,y,time)
%RIS_SEA_SURFACE Rebuild the exact incident waves used to drive RIS_float.
% time is relative to motion.parameters.startTime. No height rescaling.
    validateattributes(x,{'numeric'},{'vector','real','finite','nonempty'});
    validateattributes(y,{'numeric'},{'vector','real','finite','nonempty'});
    validateattributes(time,{'numeric'},{'scalar','real','finite','nonnegative'});
    [X,Y] = meshgrid(x,y);
    Z = reshape(RIS_wave_field(motion,[X(:),Y(:)],time),size(X));
end
