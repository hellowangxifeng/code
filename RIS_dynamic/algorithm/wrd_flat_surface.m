function sea = wrd_flat_surface()
%WRD_FLAT_SURFACE Still-water plane z = 0 in the same interface as WRD_SURFACE.
%   This reproduces the static water-air boundary used by the previous
%   "no-RIS direct" baseline, so it can be run frame by frame next to the
%   rough-sea link with the same TX motion, MC seed and receiver.
motion.parameters.startTime = 0;
motion.wave = struct('omega',0,'coefficient',0,'kx',0,'ky',0);
settings = struct('gridStep',1,'gridHalfWidth',0,'marchStep',1,'rootTolerance',1e-9);
sea = wrd_surface(motion,0,settings,[0;0]);
assert(sea.amplitudeBound==0,'wrd:FlatSurface','The flat surface must have zero amplitude.');
end
