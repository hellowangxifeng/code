function [intensity,ring] = rd_pixel_intensity(support,values,slopeX,slopeY,pixels,pixelSize,tile,neighbours)
%RD_PIXEL_INTENSITY Pixel-averaged coherent sum of tile fields (neighbour rings).
%   Tile k is known on the pixel list support{k} by its complex flux phasor
%   U_k (values{k}) and local phase slope k_k = grad arg U_k (slopeX/Y).
%   Inside a pixel of size dx x dy each tile is a local plane carrier, so
%   (RIS_phase eq. 64)
%     <|sum U|^2> = sum_k |U_k|^2
%                 + 2 sum_{k<l} Re(U_k U_l^*) sinc(dk_x dx/2pi) sinc(dk_y dy/2pi).
%   Cross terms are evaluated exactly for tiles whose array indices differ
%   by at most `neighbours` in both directions (Chebyshev ring 1..neighbours);
%   ring(r) returns the pixel-summed cross irradiance of ring r, so the
%   decay with r shows the size of the omitted rings. Negative pixel values
%   (possible only from round-off) are clipped to zero.
K = numel(support);
intensity = zeros(pixels,1);
for k = 1:K
    intensity(support{k}) = intensity(support{k})+abs(values{k}).^2;
end
ring = zeros(1,neighbours);
if K<2 || neighbours<1
    intensity = max(intensity,0);
    return
end
lookup = zeros(max(tile(:,1)),max(tile(:,2)));
lookup(sub2ind(size(lookup),tile(:,1),tile(:,2))) = 1:K;
for k = 1:K
    if isempty(support{k}), continue; end
    for dj = 0:neighbours
        for di = -neighbours:neighbours
            if dj==0 && di<=0, continue; end      % each unordered pair once
            i = tile(k,1)+di; j = tile(k,2)+dj;
            if i<1 || j<1 || i>size(lookup,1) || j>size(lookup,2), continue; end
            l = lookup(i,j);
            if l==0 || isempty(support{l}), continue; end
            [common,ia,ib] = intersect(support{k},support{l});
            if isempty(common), continue; end
            ax = (slopeX{k}(ia)-slopeX{l}(ib))*pixelSize(1)/(2*pi);
            ay = (slopeY{k}(ia)-slopeY{l}(ib))*pixelSize(2)/(2*pi);
            term = 2*real(values{k}(ia).*conj(values{l}(ib))).*normalisedSinc(ax).*normalisedSinc(ay);
            intensity(common) = intensity(common)+term;
            r = max(abs(di),dj);
            ring(r) = ring(r)+sum(term);
        end
    end
end
intensity = max(intensity,0);
end

function v = normalisedSinc(x)
v = ones(size(x));
nz = abs(x)>1e-12;
v(nz) = sin(pi*x(nz))./(pi*x(nz));
end
