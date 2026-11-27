function bank=rp_stage2_sort_bank(bank)
%RP_STAGE2_SORT_BANK Stable candidate identifiers after resumable extensions.
names={'unmodulated','steering','gaussian','flat_top','robust'};
keys=zeros(numel(bank),3);
for j=1:numel(bank)
    keys(j,:)=[find(strcmp(names,bank{j}.family),1),bank{j}.targetWidthM,bank{j}.margin];
end
[~,order]=sortrows(keys,[1 2 3]); bank=bank(order);
end
