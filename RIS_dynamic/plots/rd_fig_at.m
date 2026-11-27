function v = rd_fig_at(curve,powerDBm,P)
%RD_FIG_AT Log-linear interpolation of an ABER/OP curve at power P (dBm).
curve = double(curve(:)).';
ok = isfinite(curve);
v = NaN;
if nnz(ok)<2, return; end
v = 10.^interp1(powerDBm(ok),log10(max(curve(ok),1e-300)),P,'linear',NaN);
end
