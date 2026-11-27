function p = rd_physical(c)
%RD_PHYSICAL Fields of c that change computed numbers (checkpoint fingerprint).
% Paths, run control, plotting and post-processing statistics are excluded.
drop = {'output','run','plot','cache','root','paths','profile','metrics'};
p = rmfield(c,intersect(fieldnames(c),drop));
end
