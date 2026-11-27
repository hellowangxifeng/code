function s = rd_fig_wind_label(U)
%RD_FIG_WIND_LABEL Panel label of a sea state.
if U==0
    s = 'calm sea';
else
    s = sprintf('U = %g m/s',U);
end
end
