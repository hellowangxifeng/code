function T = rd_entry_table(E)
%RD_ENTRY_TABLE Tabulate design candidates (struct array from RD_DESIGN).
names = {'key','family','widthM','margin','valid','tiltDeg','equivalentRadiusM', ...
    'minorRadiusM','axisRatio','r50M','r90M','meanGain','aber','op'};
if isempty(E)
    T = cell2table(cell(0,numel(names)),'VariableNames',names);
    return
end
s = [E.spot];
T = table(string({E.key}).',string({E.family}).',[E.width].',[E.margin].',[E.valid].', ...
    [E.tiltDeg].',[s.equivalentRadiusM].',[s.minorRadiusM].',[s.axisRatio].', ...
    [s.r50M].',[s.r90M].',[E.meanGain].',[E.aber].',[E.op].','VariableNames',names);
end
