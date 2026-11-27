function keys = rd_frame_keys(c,m,stepS)
%RD_FRAME_KEYS Indices of the motion time samples used as optical frames.
n = round(stepS/c.ris.dt);
assert(n>=1 && abs(n*c.ris.dt-stepS)<1e-9,'rd:FrameStep', ...
    'The frame step must be a positive multiple of c.ris.dt.');
keys = 1:n:numel(m.t);
end
