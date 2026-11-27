function [tx,rx] = rd_poses(g0,n)
%RD_POSES Fixed underwater TX and nominal receiver plane for rp_geometry.
% The TX does not move: both links see exactly the same source in every frame.
tx = struct('position',repmat(g0.pTx,1,n),'direction',repmat(g0.dTx,1,n), ...
    'rotation',repmat(g0.RTx,1,1,n));
rx = struct('position',repmat(g0.rxDesired,1,n),'nominal',g0.rxDesired);
end
