function sea = rd_flat_surface()
%RD_FLAT_SURFACE Still-water plane z = 0 behind the same interface as RD_SURFACE.
% This is the calm-sea boundary of the OLD "no RIS direct" baseline; it is
% run through exactly the same Monte Carlo as the rough sea.
motion.parameters.startTime = 0;
motion.wave = struct('omega',0,'coefficient',0,'kx',0,'ky',0);
settings = struct('gridStep',1,'gridHalfWidth',0,'marchStep',1,'rootTolerance',1e-9);
sea = rd_surface(motion,0,settings,[0;0]);
assert(sea.amplitudeBound==0,'rd:FlatSurface','The flat surface must have zero amplitude.');
end
