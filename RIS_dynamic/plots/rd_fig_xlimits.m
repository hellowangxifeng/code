function x = rd_fig_xlimits(P,Y,yLimits,padding)
%RD_FIG_XLIMITS Power range in which the curves actually change inside yLimits.
%   Starts where the first curve leaves its ceiling and ends where the last
%   curve has settled on its floor (or left the bottom of the axes), padded
%   by `padding` dB and clipped to the computed power grid P.
lo = Inf; hi = -Inf;
for s = 1:size(Y,1)
    y = double(Y(s,:));
    y(~isfinite(y)) = 1;
    top = max(y);
    settled = max(1.05*y(end),yLimits(1));
    first = find(y<0.9*top & y>yLimits(1),1);
    if isempty(first), continue; end
    last = find(y>settled,1,'last');
    if isempty(last), last = first; end
    lo = min(lo,P(first));
    hi = max(hi,P(min(numel(P),last+1)));
end
if ~(lo<hi), lo = P(1); hi = P(end); end
x = [max(P(1),lo-padding),min(P(end),hi+padding)];
end
