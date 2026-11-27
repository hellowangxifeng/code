function [nodes, weights] = bsf_quadrature(order)
%BSF_QUADRATURE Gauss-Legendre nodes and weights on [-1,1].
arguments
    order (1,1) {mustBeNumeric, mustBeReal, mustBeInteger, mustBePositive}
end
j = (1:order-1).';
beta = j./sqrt(4*j.^2-1);
[vectors, values] = eig(diag(beta,1)+diag(beta,-1));
[nodes, indices] = sort(diag(values));
weights = 2*vectors(1,indices).'.^2;
end
