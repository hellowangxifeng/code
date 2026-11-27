function d = rd_design(c,z,incidence,specs,samples)
%RD_DESIGN Calm-sea RIS phase designs and receiver alignment for one height z.
%   d = rd_design(c,z,incidence,[],samples) rebuilds the RIS_phase stage-two
%   candidates on the motionless board (rp_stage2_scene + rp_phase +
%   rp_finite_phase): Steering-only and the flat_top / robust families over
%   c.study.widthsM x c.study.centerMargins. Every candidate is scored on the
%   calm sea with the SAME packet Monte Carlo and receiver (rd_capture) that
%   is used later under waves, at c.study.designPowerW. The selection rules
%   follow RIS_phase/rp_stage2_distance + rp_stage2_select:
%     * per (family, width): the margin with minimum ABER among margins whose
%       spot tilt is within c.study.directionSlackDeg of the smallest tilt;
%     * feasible: minor and equivalent radius >= c.study.minimumWidthM,
%       axis ratio <= c.study.maximumAxisRatio, r50 >= 2 cm, tilt within slack;
%     * Robust ABER / Robust OP: best feasible robust-family design;
%       Max mean power: best feasible design of either family;
%       Flat top: flat_top design closest to c.study.referenceWidthM.
%   d = rd_design(c,z,incidence,specs,samples) skips the search and builds the
%   listed specs (struct array with key, family, width, margin), as used by
%   the incidence sweep (nominal widths kept, carriers recomputed, as Fig10).
%
%   Each receiver is aligned once, at the calm-sea ballistic centroid of its
%   own link (RIS_phase convention). If the calm direct link carries no
%   ballistic power (total internal reflection), d.direct.center is NaN and
%   rd_run aligns it with a rough-sea pilot run (rd_pilot).
fprintf('[design] z = %.2f m, incidence = %.2f deg\n',z,incidence);
[g,~,a,cc] = rp_stage2_scene(c,z,incidence);
Z = g.rxPlaneZ;
baseline = rp_phase(cc,g,a);
transport = rd_transport(cc,cc.dynamic.photons);
[~,packets] = wris_mc(transport,g,cc.mc.seed);
lcp = rd_ris_packets(cc,g,packets);
air = rd_direct_mc(transport,g,rd_flat_surface(),cc.mc.seed);
top = g.pRis+g.RRis*[0;0;cc.ris.thickness];
s = -g.pTx(3)/g.dTx(3);
crossing = g.pTx+s*g.dTx;
d = struct('distanceM',z,'incidenceDeg',incidence,'Z',Z,'g',g, ...
    'originRis',top.','originDirect',[crossing(1:2).',0],'crossXY',crossing(1:2).');
spot = rd_spot(rd_to_plane(air,Z));
center = spot.centerXY;
alignment = 'calm';
if ~(spot.ballisticPowerW>=cc.receiver.minimumCalmPowerW)
    center = [NaN NaN];
    alignment = 'pilot required';
end
d.direct = struct('spot',spot,'center',center,'alignment',alignment, ...
    'transmittedPowerW',air.transmittedPowerW,'tirFraction',air.tirFraction);

if isempty(specs)
    list = struct('key','steering','family','steering','width',0,'margin',0);
    for family = {'flat_top','robust'}
        for width = cc.study.widthsM
            for margin = cc.study.centerMargins
                list(end+1) = struct('key','','family',family{1},'width',width, ...
                    'margin',margin); %#ok<AGROW>
            end
        end
    end
else
    list = specs;
end
E = repmat(emptyEntry(),1,numel(list));
for j = 1:numel(list)
    E(j) = candidate(cc,g,a,baseline,lcp,Z,d.originRis,samples,list(j));
    fprintf('  %3d/%3d %-9s w = %5.3f m, margin = %5.3f | R_eq = %6.2f cm | ABER = %.3g, OP = %.3g %s\n', ...
        j,numel(list),E(j).family,E(j).width,E(j).margin,100*E(j).spot.equivalentRadiusM, ...
        E(j).aber,E(j).op,E(j).message);
end
d.screening = rd_entry_table(E);

if ~isempty(specs)
    d.ris = E;
    d.bank = E([]);
    d.bankFeasible = false(1,0);
    d.bankTable = rd_entry_table(d.bank);
    return
end

% ----- RIS_phase stage-two selection on the calm sea -----
family = string({E.family});
width = [E.width];
valid = [E.valid];
bank = zeros(1,0);
for f = ["flat_top","robust"]
    for w = cc.study.widthsM
        ids = find(family==f & abs(width-w)<1e-12 & valid);
        if isempty(ids), continue; end
        tilt = [E(ids).tiltDeg];
        ids = ids(tilt<=min(tilt)+cc.study.directionSlackDeg);
        [~,best] = min([E(ids).aber]);
        bank(end+1) = ids(best); %#ok<AGROW>
    end
end
assert(~isempty(bank),'rd:Design','No valid flat_top/robust candidate at z = %g m.',z);
B = E(bank);
spots = [B.spot];
feasible = [spots.minorRadiusM]>=cc.study.minimumWidthM & ...
    [spots.equivalentRadiusM]>=cc.study.minimumWidthM & ...
    [spots.axisRatio]<=cc.study.maximumAxisRatio & [spots.r50M]>=0.02;
if any(feasible)
    tilt = [B.tiltDeg];
    feasible = feasible & tilt<=min(tilt(feasible))+cc.study.directionSlackDeg;
else
    warning('rd:Feasible','z = %g m: no candidate meets the RIS_phase spot constraints; using all valid ones.',z);
    feasible = true(size(B));
end
bankFamily = string({B.family});
robust = feasible & bankFamily=="robust";
if ~any(robust)
    warning('rd:Robust','z = %g m: no feasible robust design; using all robust designs.',z);
    robust = bankFamily=="robust";
end
score = [B.aber]; score(~robust) = Inf; [~,iAber] = min(score);
score = [B.op]; score(~robust) = Inf; [~,iOp] = min(score);
score = [B.meanGain]; score(~feasible) = -Inf; [~,iPower] = max(score);
flat = find(bankFamily=="flat_top");
if isempty(flat)
    warning('rd:FlatTop','z = %g m: no valid flat-top design; Flat top uses Robust ABER.',z);
    iFlat = iAber;
else
    [~,nearest] = min(abs([B(flat).width]-cc.study.referenceWidthM));
    iFlat = flat(nearest);
end
steering = find(family=="steering",1);
picked = [E(steering),B(iFlat),B(iAber),B(iOp),B(iPower)];
keys = {'steering','flat_top','robust_aber','robust_op','max_power'};
for k = 1:numel(keys), picked(k).key = keys{k}; end
d.ris = picked;
d.bank = B;
d.bankFeasible = feasible;
d.bankTable = rd_entry_table(B);
d.bankTable.feasible = feasible(:);
end

function e = candidate(c,g,a,baseline,lcp,Z,origin,samples,spec)
e = emptyEntry();
e.key = spec.key; e.family = spec.family; e.width = spec.width; e.margin = spec.margin;
try
    p = rp_finite_phase(c,g,a,spec.family,spec.width,spec.margin,baseline);
catch err
    e.message = sprintf('(%s)',err.message);
    return
end
e.p = rd_slim_phase(p,spec.family,spec.width,spec.margin);
land = rd_trace(c,g,lcp,e.p,Z);
e.spot = rd_spot(land);
if ~(e.spot.ballisticPowerW>0)
    e.message = '(no ballistic power on the receiver plane)';
    return
end
e.center = e.spot.centerXY;
e.tiltDeg = atan2d(norm(e.center-g.pRis(1:2).'),Z-g.pRis(3)-c.ris.thickness);
rx = rd_receiver(c,origin,e.center,Z,samples);
gain = rd_capture(c,rd_compress(land,c.compress),rx);
m = rd_metrics(c,gain,rx,c.study.designPowerW);
e.meanGain = m.meanGain; e.aber = m.aber; e.op = m.op; e.valid = true;
end

function e = emptyEntry()
e = struct('key','','family','','width',NaN,'margin',NaN,'valid',false,'p',[], ...
    'spot',rd_spot([]),'center',[NaN NaN],'tiltDeg',NaN,'meanGain',NaN,'aber',NaN, ...
    'op',NaN,'message','');
end
