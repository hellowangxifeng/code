function links = wrd_fig_links(results)
%WRD_FIG_LINKS Display order, labels and colours of the links present in results.
%   key matches the frame field: 'flat' (still-water baseline, optional),
%   'direct' (rough-sea direct MC) and 'ris' (moving RIS).
links = struct('key',{'flat','direct','ris'}, ...
    'label',{'Still-water direct (old baseline)','Rough-sea direct (MC)','Moving RIS'}, ...
    'color',{[0.45 0.45 0.45],[0.74 0.29 0.19],[0.10 0.36 0.55]}, ...
    'style',{'--','-','-'},'marker',{'none','.','.'});
f = results.cases{1}.frames;
if ~isfield(f,'flat') || isempty(f(1).flat)
    links(1) = [];
end
end
