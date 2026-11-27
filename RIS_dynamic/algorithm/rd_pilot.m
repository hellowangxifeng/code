function centers = rd_pilot(c,wind,g0,crossXY,Z)
%RD_PILOT Align a direct-link receiver that has no calm-sea ballistic light.
% Beyond the critical angle the calm surface transmits no ballistic power,
% so the calm centroid does not exist. The receiver is then aligned at the
% long-term (time-averaged) ballistic centroid of the rough-sea direct link,
% estimated from c.receiver.pilotFrames frames with independent MC seeds.
fprintf('[pilot] direct-link receiver alignment on the rough sea, U = %g m/s\n',wind);
m = rd_motion(c,wind,c.ris.seed);
keys = unique(round(linspace(1,numel(m.t),c.receiver.pilotFrames)));
transport = rd_transport(c,c.receiver.pilotPhotons);
[tx,rxPose] = rd_poses(g0,numel(m.t));
nz = numel(Z);
ballistic = zeros(nz,3);
total = zeros(nz,3);
for f = 1:numel(keys)
    g = rp_geometry(c,m,tx,rxPose,keys(f));
    sea = rd_surface(m,m.t(keys(f)),c.surface,crossXY);
    air = rd_direct_mc(transport,g,sea,c.mc.seed+900000+f);
    for iz = 1:nz
        land = rd_to_plane(air,Z(iz));
        b = land.ballistic;
        ballistic(iz,:) = ballistic(iz,:)+[sum(land.power(b)),land.power(b).'*land.xy(b,:)];
        total(iz,:) = total(iz,:)+[sum(land.power),land.power.'*land.xy];
    end
end
centers = nan(nz,2);
for iz = 1:nz
    if ballistic(iz,1)>0
        centers(iz,:) = ballistic(iz,2:3)/ballistic(iz,1);
    elseif total(iz,1)>0
        centers(iz,:) = total(iz,2:3)/total(iz,1);
    end
end
end
