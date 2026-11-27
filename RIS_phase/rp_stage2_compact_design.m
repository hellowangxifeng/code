function design=rp_stage2_compact_design(d)
%RP_STAGE2_COMPACT_DESIGN Send only phase descriptors to incidence workers.
bank=cell(1,6);
for k=3:6
    b=d.bank{d.comparisonIndices(k)};
    bank{k}=struct('family',b.family,'targetWidthM',b.targetWidthM,'margin',b.margin);
end
design=struct('distanceM',d.distanceM,'bank',{bank},'comparisonIndices',1:6, ...
    'selectedCandidates',d.comparisonIndices);
end
