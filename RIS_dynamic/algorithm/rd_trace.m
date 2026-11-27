function pk = rd_trace(c,g,L,p,Z)
%RD_TRACE LC packets -> frozen tile phase ramps -> glass/air -> plane z = Z.
%   Uses RIS_phase/rp_trace unchanged with the CURRENT board pose g and the
%   phase p designed on the motionless board (rd_slim_phase).
g.rxPlaneZ = Z;
lookup = zeros(max([p.unitIds(:);L.unitId(:);1]),1);
lookup(p.unitIds) = 1:numel(p.unitIds);
m = lookup(L.unitId);
keep = m>0;
a = struct('lcXY',L.lcXY(keep,:),'qIn',L.qIn(keep,:),'powerW',L.powerW(keep), ...
    'unitId',L.unitId(keep),'sampleUnitIndex',double(m(keep)));
r = rp_trace(c,g,a,p);
ok = r.valid & r.powerW>0;
ballistic = L.ballistic(keep);
pk = struct('xy',r.landing(ok,1:2),'dir',r.directionWorld(ok,:),'power',r.powerW(ok), ...
    'ballistic',ballistic(ok));
end
