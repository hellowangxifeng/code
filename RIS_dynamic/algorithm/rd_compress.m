function pk = rd_compress(pk,settings)
%RD_COMPRESS Merge packets that share a fine position cell and direction cell.
% Power is conserved exactly; merged positions/directions are power-weighted
% means, so the quantisation error is below settings.positionM/2 and
% settings.directionCos/2 (1 mm / 1 mrad for the defaults), far below the
% 5 cm aperture radius and the 130 mrad FOV.
if numel(pk.power)<2, return; end
key = [round(pk.xy/settings.positionM),round(pk.dir(:,1:2)/settings.directionCos), ...
    double(pk.ballistic)];
[~,~,group] = unique(key,'rows');
w = accumarray(group,pk.power);
xy = [accumarray(group,pk.power.*pk.xy(:,1)),accumarray(group,pk.power.*pk.xy(:,2))]./w;
d = [accumarray(group,pk.power.*pk.dir(:,1)),accumarray(group,pk.power.*pk.dir(:,2)), ...
    accumarray(group,pk.power.*pk.dir(:,3))]./w;
d = d./vecnorm(d,2,2);
b = accumarray(group,double(pk.ballistic),[],@max)>0;
pk = struct('xy',xy,'dir',d,'power',w,'ballistic',b);
end
