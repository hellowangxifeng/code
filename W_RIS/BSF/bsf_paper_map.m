function result = bsf_paper_map(c,g,x,y,order,frequencyNodes)
%BSF_PAPER_MAP Cell-averaged (3-18) with matched finite Gaussian waist.
arguments
    c (1,1) struct
    g (1,1) struct
    x {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
    y {mustBeNumeric,mustBeReal,mustBeFinite,mustBeVector}
    order (1,1) {mustBeInteger,mustBePositive,mustBeReal} = 3
    frequencyNodes (1,1) {mustBeInteger,mustBePositive,mustBeReal} = 2049
end
[X,Y] = ndgrid(x,y);
[nodes,weights] = bsf_quadrature(order);
dx = mean(diff(x));
dy = mean(diff(y));
Ib = zeros(size(X));
Is = Ib;
for i = 1:order
    for j = 1:order
        point = bsf_paper(c,g,X+dx/2*nodes(i),Y+dy/2*nodes(j),FrequencyNodes=frequencyNodes);
        Ib = Ib+weights(i)*weights(j)/4*point.Ib;
        Is = Is+weights(i)*weights(j)/4*point.Is;
    end
end
result = struct('Ib',Ib,'Is',Is,'IBSF',Ib+Is,'areaOrder',order, ...
    'frequencyNodes',frequencyNodes,'minimumScatter',min(Is,[],'all'));
end
