function idx = rd_fig_index(values,wanted)
%RD_FIG_INDEX Index of the stored value closest to each wanted value.
idx = zeros(size(wanted));
for j = 1:numel(wanted)
    [~,idx(j)] = min(abs(values(:)-wanted(j)));
end
end
