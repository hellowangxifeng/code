function [xx,yy] = atlas_pchip(x,y,n)
%ATLAS_PCHIP Shape-preserving dense interpolation used for all line plots.
if nargin<3, n=300; end
x = x(:); y = y(:);
[x,idx] = sort(x); y = y(idx);
[xu,ia] = unique(x,'stable'); yu = y(ia);
if numel(xu) < 2
    xx = xu; yy = yu; return;
end
xx = linspace(xu(1),xu(end),n).';
yy = interp1(xu,yu,xx,'pchip');
end
