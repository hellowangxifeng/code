function s = rd_fig_style()
%RD_FIG_STYLE Colours / markers of the five links (Direct, Steering, Flat, Robust, Max power).
s.colors = [.40 .40 .40; .05 .35 .72; .90 .58 .04; .10 .56 .25; .83 .16 .16];
s.lines = {':','-','-','-','-'};
s.markers = {'none','s','^','o','d'};
end
