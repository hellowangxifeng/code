function p = rd_physical(c)
%RD_PHYSICAL Fields of c that change computed numbers (checkpoint fingerprint).
% Paths, run control, plotting and post-processing statistics are excluded.
% The wind LIST is excluded too: every condition, training and pilot
% fingerprint records its own wind, so adding or removing a wind speed
% never invalidates the calm-sea designs or the other winds.
drop = {'output','run','plot','cache','root','paths','profile','metrics'};
p = rmfield(c,intersect(fieldnames(c),drop));
if isfield(p,'sea') && isfield(p.sea,'windsMps')
    p.sea = rmfield(p.sea,'windsMps');
end
end
