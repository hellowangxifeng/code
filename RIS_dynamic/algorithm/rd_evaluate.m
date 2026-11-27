function R = rd_evaluate(c,designs,training,main,incidence)
%RD_EVALUATE Turn stored instantaneous gains into ABER / OP / CDF / tables.
%   Light result structure R (read by the plots; no packets, no phases):
%     R.main.aber, R.main.op       [wind x z x link x turbulence x power]
%     R.main.meanGain, zeroFraction, scintillation, medianGainDB  [wind x z x link]
%     R.main.cdf                   [wind x z x link x edge], gain in dB
%     R.main.frames{wind}          time series (tilt, slope, TIR, frame ABER/OP)
%     R.main.maps{wind}            receiver-plane irradiance (calm / snapshot / mean)
%     R.incidence.aber, .op        [angle x {calm, wind} x z x link x power]
%     R.summary, R.designTable, R.trainingTable, R.incidenceTable
%   Links follow rd_schemes(); winds include the calm sea (U = 0) first.
sch = rd_schemes();
if strcmpi(c.training.objective,'op'), sch.names{end} = 'Sea-aware OP'; end
if ~c.training.enabled, sch.names{end} = 'Sea-aware (= Robust ABER)'; end
K = numel(sch.keys);
P = c.study.powerDBm(:).';
powerW = 10.^(P/10)/1000;
np = numel(P);
turbulence = c.metrics.turbulence;
nt = size(turbulence,1);
z = c.scene.distancesM;
nz = numel(z);
winds = cellfun(@(b) b.wind,main);
nw = numel(winds);
edges = c.metrics.cdfEdgesDB(:).';
fixed = c.metrics.fixedPowerDBm;

R = struct('config',c,'schemes',sch,'powerDBm',P,'distancesM',z,'winds',winds, ...
    'turbulence',{turbulence});
M = struct();
M.aber = nan(nw,nz,K,nt,np); M.op = M.aber;
M.meanGain = nan(nw,nz,K); M.zeroFraction = M.meanGain; M.scintillation = M.meanGain;
M.medianGainDB = M.meanGain; M.centroidWanderM = M.meanGain;
M.cdfEdgesDB = edges; M.cdf = nan(nw,nz,K,numel(edges));
M.fovFloor = nan(nz,K); M.receiverCenter = nan(nz,K,2);
M.frames = cell(1,nw); M.maps = cell(1,nw);
M.Hs = nan(1,nw); M.maxTiltDeg = nan(1,nw); M.smallAngleValid = true(1,nw);
rows = {};
for iw = 1:nw
    b = main{iw};
    nf = numel(b.frames.t);
    fr = b.frames;
    fr.meanGainFrames = nan(nf,nz,K);
    fr.aberFrames = nan(nf,np,nz,K,'single');
    fr.opFrames = fr.aberFrames;
    fr.planePowerW = b.planePowerW;
    fr.centroidXY = b.centroidXY;
    for iz = 1:nz
        for k = 1:K
            rx = b.receivers{iz,k};
            if iw==1
                M.fovFloor(iz,k) = rx.fovOutageExact;
                M.receiverCenter(iz,k,:) = reshape(rx.center,1,1,2);
            end
            G = b.gain(:,:,iz,k);
            if ~rx.available || all(isnan(G(:))), continue; end
            for it = 1:nt
                cc = c;
                cc.study.turbulence = turbulence{it,1};
                cc.study.atmosphere = turbulence{it,2};
                mm = rd_metrics(cc,G,rx,powerW);
                M.aber(iw,iz,k,it,:) = mm.aber;
                M.op(iw,iz,k,it,:) = mm.op;
                if it==1
                    fr.aberFrames(:,:,iz,k) = single(mm.aberFrames);
                    fr.opFrames(:,:,iz,k) = single(mm.opFrames);
                    fr.meanGainFrames(:,iz,k) = mm.meanGainFrames;
                end
                rows(end+1,:) = {winds(iw),b.Hs,z(iz),string(sch.keys{k}),string(sch.names{k}), ...
                    it,string(turbulence{it,1}),string(turbulence{it,2}), ...
                    10*log10(max(mm.meanGain,realmin)),NaN,NaN,NaN, ...
                    at(mm.aber,P,fixed),at(mm.op,P,fixed)}; %#ok<AGROW>
            end
            g = double(G(:));
            dB = 10*log10(g);
            M.meanGain(iw,iz,k) = mean(g);
            M.zeroFraction(iw,iz,k) = mean(g<=0);
            M.scintillation(iw,iz,k) = var(g)/max(mean(g)^2,realmin);
            M.medianGainDB(iw,iz,k) = median(dB);
            M.cdf(iw,iz,k,:) = reshape(mean(dB<=edges,1),1,1,1,[]);
            xy = reshape(b.centroidXY(:,:,iz,k),nf,2);
            xy = xy(all(isfinite(xy),2),:);
            if size(xy,1)>1, M.centroidWanderM(iw,iz,k) = sqrt(sum(var(xy,1,1))); end
            for r = size(rows,1)-nt+1:size(rows,1)
                rows{r,10} = M.medianGainDB(iw,iz,k);
                rows{r,11} = M.zeroFraction(iw,iz,k);
                rows{r,12} = M.scintillation(iw,iz,k);
            end
        end
    end
    M.frames{iw} = fr;
    M.maps{iw} = struct('mean',b.meanMap,'snapshot',b.snapshotMap,'snapshotTimeS',b.snapshotTimeS, ...
        'halfWidthM',b.mapHalfWidthM,'bins',c.maps.bins);
    M.Hs(iw) = b.Hs;
    if isfield(b.validity,'maxTiltDeg'), M.maxTiltDeg(iw) = b.validity.maxTiltDeg; end
    if isfield(b.validity,'smallAngleUnder10Deg'), M.smallAngleValid(iw) = b.validity.smallAngleUnder10Deg; end
end
R.main = M;
names = {'windMps','HsM','distanceM','scheme','name','turbulenceId','ocean','atmosphere', ...
    'meanGainDB','medianGainDB','zeroFraction','scintillation','aberFixed','opFixed'};
if isempty(rows), rows = cell(0,numel(names)); end
R.summary = cell2table(rows,'VariableNames',names);
R.summary.Properties.Description = sprintf('ABER/OP at %g dBm; gains pooled over frames and jitter',fixed);

% ----- Designs -----
designRows = {};
for iz = 1:nz
    d = designs{iz};
    T = rd_entry_table(d.ris);
    T.distanceM = repmat(z(iz),height(T),1);
    designRows{end+1} = T; %#ok<AGROW>
end
R.designTable = vertcat(designRows{:});
R.designs = cellfun(@(d) struct('distanceM',d.distanceM,'Z',d.Z,'screening',d.screening, ...
    'bankTable',d.bankTable,'selected',rd_entry_table(d.ris),'direct',d.direct), ...
    designs,'UniformOutput',false);

% ----- Sea-aware training -----
trainingRows = {};
for iw = 1:numel(training)
    T = training{iw};
    for iz = 1:numel(T)
        for j = 1:numel(T(iz).bankIds)
            trainingRows(end+1,:) = {T(iz).windMps,T(iz).distanceM,string(T(iz).family{j}), ...
                T(iz).widthM(j),T(iz).margin(j),T(iz).equivalentRadiusM(j),T(iz).aber(j), ...
                T(iz).op(j),T(iz).meanGain(j),T(iz).bankIds(j)==T(iz).chosen}; %#ok<AGROW>
        end
    end
end
if isempty(trainingRows), trainingRows = cell(0,10); end
R.trainingTable = cell2table(trainingRows,'VariableNames',{'windMps','distanceM','family', ...
    'widthM','margin','equivalentRadiusM','aber','op','meanGain','chosen'});

% ----- Incidence sweep -----
R.incidence = [];
R.incidenceTable = table();
if ~isempty(incidence.anglesDeg)
    A = incidence.anglesDeg;
    na = numel(A);
    I = struct('anglesDeg',A,'windMps',incidence.windMps,'aber',nan(na,2,nz,K,np), ...
        'op',nan(na,2,nz,K,np),'meanGain',nan(na,2,nz,K),'zeroFraction',nan(na,2,nz,K), ...
        'ballisticTirFraction',nan(na,2),'directTransmittedW',nan(na,2), ...
        'criticalAngleDeg',asind(c.optics.nAir/c.optics.nWater));
    incRows = {};
    for ia = 1:na
        for k2 = 1:2
            b = incidence.blocks{ia,k2};
            I.ballisticTirFraction(ia,k2) = mean(b.frames.ballisticTirFraction);
            I.directTransmittedW(ia,k2) = mean(b.frames.directTransmittedW);
            for iz = 1:nz
                for k = 1:K
                    rx = b.receivers{iz,k};
                    G = b.gain(:,:,iz,k);
                    if ~rx.available || all(isnan(G(:))), continue; end
                    mm = rd_metrics(c,G,rx,powerW);
                    I.aber(ia,k2,iz,k,:) = mm.aber;
                    I.op(ia,k2,iz,k,:) = mm.op;
                    I.meanGain(ia,k2,iz,k) = mm.meanGain;
                    I.zeroFraction(ia,k2,iz,k) = mean(G(:)<=0);
                    incRows(end+1,:) = {A(ia),b.wind,z(iz),string(sch.keys{k}), ...
                        10*log10(max(mm.meanGain,realmin)),I.zeroFraction(ia,k2,iz,k), ...
                        at(mm.aber,P,fixed),at(mm.op,P,fixed)}; %#ok<AGROW>
                end
            end
        end
    end
    R.incidence = I;
    if ~isempty(incRows)
        R.incidenceTable = cell2table(incRows,'VariableNames',{'incidenceDeg','windMps', ...
            'distanceM','scheme','meanGainDB','zeroFraction','aberFixed','opFixed'});
    end
    R.incidenceDesigns = incidence.designs;
end
end

function v = at(curve,P,fixed)
curve = double(curve(:)).';
ok = isfinite(curve);
v = NaN;
if nnz(ok)<2, return; end
v = 10.^interp1(P(ok),log10(max(curve(ok),1e-300)),fixed,'linear',NaN);
end
