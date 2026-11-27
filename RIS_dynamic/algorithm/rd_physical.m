function p = rd_physical(c)
%RD_PHYSICAL Fields of c that change computed numbers (checkpoint fingerprint).
% Paths, run control, plotting, post-processing statistics and the grid
% LISTS are excluded: every condition / design fingerprint records its own
% wind and height, so adding a wind or a height never invalidates the others.
drop = {'output','run','plot','cache','root','paths','profile','metrics','rank','main'};
p = rmfield(c,intersect(fieldnames(c),drop));
p.sea = rmfield(p.sea,intersect(fieldnames(p.sea),{'windsMps'}));
p.scene = rmfield(p.scene,intersect(fieldnames(p.scene),{'distancesM'}));
p.study = rmfield(p.study,intersect(fieldnames(p.study),{'turbulence','atmosphere','powerDBm'}));
end
