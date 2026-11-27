function s = rd_schemes()
%RD_SCHEMES Link order used in every result array (dimension "k").
s.keys = {'direct','steering','flat_top','robust_aber','robust_op','max_power','sea_aware'};
s.names = {'No RIS direct','Steering-only','Flat top','Robust ABER','Robust OP', ...
    'Max mean power','Sea-aware ABER'};
end
