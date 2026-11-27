function [floorW,topW] = wrd_fig_range(results,links)
%WRD_FIG_RANGE Common log power axis for every tile and link (W).
values = [];
tops = [];
for id = 1:numel(results.cases)
    f = results.cases{id}.frames;
    for j = 1:numel(links)
        v = [f.(links(j).key)];
        values = [values,[v.rxPowerW]]; %#ok<AGROW>
        tops = [tops,[v.rxPowerW]+1.96*[v.powerSEW]]; %#ok<AGROW>
    end
end
positive = values(values>0);
floorW = results.config.metrics.displayFloorW;
if ~isempty(positive), floorW = min(positive)/3; end
floorW = min(floorW,results.config.metrics.outageThresholdW/3);
topW = 1.25*max([floorW*50,tops,results.config.metrics.outageThresholdW*3]);
end
