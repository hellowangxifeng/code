function m = rd_static_motion(g0,n)
%RD_STATIC_MOTION Motionless board at its flotation equilibrium (calm sea).
m = struct('t',zeros(1,n),'rotationMatrix',repmat(g0.RRis,1,1,n), ...
    'bottomCenter',repmat(g0.pRis,1,n));
end
