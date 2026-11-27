function label = wrd_fig_title(item)
%WRD_FIG_TITLE "U = 5.0 m/s, theta_0 = 40.0 deg" for one case.
if isfield(item,'nominalIncidenceDeg')
    incidence = item.nominalIncidenceDeg;
else
    a = item.txRotationXYDeg;
    incidence = acosd(cosd(a(1))*cosd(a(2)));
end
label = sprintf('U = %.1f m/s, \\theta_0 = %.1f^{\\circ}',item.windMps,incidence);
end
