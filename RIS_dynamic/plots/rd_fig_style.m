function s = rd_fig_style()
%RD_FIG_STYLE Colours and line styles per link, in rd_schemes order.
% The first six colours are RIS_phase/rp_s2_style; sea-aware is new.
s.colors = [.35 .35 .35; .05 .35 .72; .90 .58 .04; .10 .56 .25; .55 .22 .67; ...
    .83 .16 .16; .00 .62 .65];
s.lines = {'-','-','-','-','--',':','-.'};
s.markers = {'o','s','^','v','d','x','p'};
s.widths = [1.6 1.2 1.2 1.2 1.2 1.4 1.4];
end
