function pk = rd_trace(c,g,L,p,Z)
%RD_TRACE MC packets on the LC screen -> loaded tile phases -> glass/air -> plane z = Z.
%   Geometric optics of RIS_phase/rp_trace (local phase gradient, Snell and
%   Fresnel at the glass/air face) for the frozen pose g and the design p.
%   L comes from RD_RIS_PACKETS (raw, unweighted MC power); ballistic and
%   scattered packets are both traced, so the result is the full-MC
%   counterpart of the complex field of RD_WAVE.
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
