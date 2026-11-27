function pk = rd_to_plane(air,Z)
%RD_TO_PLANE Straight air propagation of surface-exit packets to the plane z = Z.
dz = air.dir(:,3);
s = (Z-air.pos(:,3))./dz;
ok = dz>1e-9 & s>0 & air.power>0;
pk = struct('xy',air.pos(ok,1:2)+air.dir(ok,1:2).*s(ok),'dir',air.dir(ok,:), ...
    'power',air.power(ok),'ballistic',air.ballistic(ok));
end
