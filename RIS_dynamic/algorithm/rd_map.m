function I = rd_map(pk,xEdges,yEdges)
%RD_MAP Horizontal irradiance (W/m^2) of landed packets on a rectangular grid.
% Rows = y cells, columns = x cells; packets outside the edges are ignored.
nx = numel(xEdges)-1;
ny = numel(yEdges)-1;
I = zeros(ny,nx);
if isempty(pk.power), return; end
ix = discretize(pk.xy(:,1),xEdges);
iy = discretize(pk.xy(:,2),yEdges);
in = ~isnan(ix) & ~isnan(iy);
if any(in)
    I = accumarray([iy(in),ix(in)],pk.power(in),[ny nx]) / ...
        ((xEdges(2)-xEdges(1))*(yEdges(2)-yEdges(1)));
end
end
