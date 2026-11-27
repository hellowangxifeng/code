function [nodes,weights] = rp_gauss_legendre(order)
%RP_GAUSS_LEGENDRE Gaussian quadrature on [-1,1].
arguments
    order (1,1) double {mustBeInteger,mustBePositive}
end
j=(1:order-1).';
beta=j./sqrt(4*j.^2-1);
[V,D]=eig(diag(beta,1)+diag(beta,-1));
[nodes,indices]=sort(diag(D));
weights=2*V(1,indices).'.^2;
end
