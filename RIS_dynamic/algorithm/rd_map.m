function I = rd_map(pk,center,halfWidth,bins)
%RD_MAP Horizontal irradiance (W/m^2) on a square grid centred on the receiver.
% Rows = y, columns = x, cell centres at ((1:bins)-0.5)/bins*2*halfWidth-halfWidth.
I = zeros(bins,bins);
if isempty(pk.power) || any(~isfinite(center)), return; end
edges = linspace(-halfWidth,halfWidth,bins+1);
ix = discretize(pk.xy(:,1)-center(1),edges);
iy = discretize(pk.xy(:,2)-center(2),edges);
in = ~isnan(ix) & ~isnan(iy);
if any(in)
    I = accumarray([iy(in),ix(in)],pk.power(in),[bins bins])/(edges(2)-edges(1))^2;
end
end
