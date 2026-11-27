function T = rd_train(c,wind,designs,samples)
%RD_TRAIN Sea-aware width selection on an INDEPENDENT sea realization.
%   The calm-sea design bank (one margin-selected candidate per family and
%   width, RD_DESIGN) is re-scored under the moving board of a sea generated
%   with wave seed c.ris.seed + c.training.seedOffset, i.e. different wave
%   phases from the evaluation sea. The candidate with the smallest
%   time-averaged ABER (or OP, c.training.objective) at c.study.designPowerW
%   becomes the "sea-aware" link for this wind. Ties go to the larger mean
%   gain. Only the RIS link is needed here.
fprintf('[training] sea-aware selection, U = %g m/s\n',wind);
m = rd_motion(c,wind,c.ris.seed+c.training.seedOffset);
keys = rd_frame_keys(c,m,c.training.frameStepS);
transport = rd_transport(c,c.training.photons);
[tx,rxPose] = rd_poses(designs{1}.g,numel(m.t));
nz = numel(designs);
ids = cell(1,nz); rx = cell(1,nz); G = cell(1,nz);
for iz = 1:nz
    d = designs{iz};
    use = ismember({d.bank.family},c.training.families) & [d.bank.valid];
    if any(use & d.bankFeasible), use = use & d.bankFeasible; end
    ids{iz} = find(use);
    rx{iz} = cell(1,numel(ids{iz}));
    for j = 1:numel(ids{iz})
        rx{iz}{j} = rd_receiver(c,d.originRis,d.bank(ids{iz}(j)).center,d.Z,samples);
    end
    G{iz} = zeros(samples.n,numel(keys),numel(ids{iz}),'single');
end
clock = tic;
for f = 1:numel(keys)
    g = rp_geometry(c,m,tx,rxPose,keys(f));
    [~,packets] = wris_mc(transport,g,c.mc.seed+700000+f);
    lcp = rd_ris_packets(c,g,packets);
    for iz = 1:nz
        d = designs{iz};
        for j = 1:numel(ids{iz})
            land = rd_trace(c,g,lcp,d.bank(ids{iz}(j)).p,d.Z);
            G{iz}(:,f,j) = single(rd_capture(c,rd_compress(land,c.compress),rx{iz}{j}));
        end
    end
    if f==1 || f==numel(keys) || mod(f,c.run.progressEvery)==0
        fprintf('    [training] U = %g m/s | frame %d/%d | %.1f min\n',wind,f,numel(keys),toc(clock)/60);
    end
end
T = repmat(struct('windMps',wind,'distanceM',NaN,'bankIds',[],'family',{{}},'widthM',[], ...
    'margin',[],'equivalentRadiusM',[],'aber',[],'op',[],'meanGain',[],'chosen',NaN, ...
    'objective',lower(c.training.objective)),1,nz);
for iz = 1:nz
    d = designs{iz};
    n = numel(ids{iz});
    aber = nan(1,n); op = aber; meanGain = aber;
    for j = 1:n
        mm = rd_metrics(c,G{iz}(:,:,j),rx{iz}{j},c.study.designPowerW);
        aber(j) = mm.aber; op(j) = mm.op; meanGain(j) = mm.meanGain;
    end
    b = d.bank(ids{iz});
    T(iz).distanceM = d.distanceM;
    T(iz).bankIds = ids{iz};
    T(iz).family = {b.family};
    T(iz).widthM = [b.width];
    T(iz).margin = [b.margin];
    if n>0
        spots = [b.spot];
        T(iz).equivalentRadiusM = [spots.equivalentRadiusM];
    end
    T(iz).aber = aber; T(iz).op = op; T(iz).meanGain = meanGain;
    if n>0 && any(isfinite(aber))
        score = aber;
        if strcmpi(c.training.objective,'op'), score = op; end
        score(~isfinite(score)) = Inf;
        best = find(score==min(score));
        [~,j] = max(meanGain(best));
        T(iz).chosen = ids{iz}(best(j));
        fprintf('  U = %g m/s, z = %g m: sea-aware = %s, w = %.3f m (ABER %.3g, OP %.3g)\n', ...
            wind,d.distanceM,b(best(j)).family,b(best(j)).width,aber(best(j)),op(best(j)));
    end
end
end
