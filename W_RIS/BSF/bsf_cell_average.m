function values = bsf_cell_average(fun, x, y, order)
%BSF_CELL_AVERAGE Area averages, with matrix convention I(ix,iy) of wris_mc.
arguments
    fun (1,1) function_handle
    x {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
    y {mustBeNumeric, mustBeReal, mustBeFinite, mustBeVector}
    order (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive} = 4
end
assert(numel(x)>1 && numel(y)>1, 'bsf:Grid', 'At least two cell centres are required.');
dx = mean(diff(x));
dy = mean(diff(y));
assert(dx>0 && dy>0 && max(abs(diff(x)-dx))<1e-10 && max(abs(diff(y)-dy))<1e-10, ...
    'bsf:Grid', 'The MC comparison requires a uniform increasing grid.');
[X,Y] = ndgrid(x,y);
[nodes,weights] = bsf_quadrature(order);
values = zeros(size(X));
for i = 1:order
    for j = 1:order
        values = values + weights(i)*weights(j)/4*fun(X+dx/2*nodes(i),Y+dy/2*nodes(j));
    end
end
end
