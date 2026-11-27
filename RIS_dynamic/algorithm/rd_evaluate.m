function R = rd_evaluate(c,C)
%RD_EVALUATE Channel statistics, scheme selection, rankings and the light result R.
%   C{iz,iw} from rd_condition. For every condition and every fading row
%   {ocean, atmosphere} of c.metrics.turbulence:
%     * ABER, OP, <P_r> of every design from its complex-field gains and of
%       the direct link from its MC gains (rd_link_metrics, same samples);
%     * the design used by each scheme at each power (rd_select);
%     * curves of the five schemes and the required-power ranking (rd_rank).
%   For the main fading row the same selections are also evaluated with the
%   MC gains of the RIS designs (validation curves).
P = c.metrics.powerDBm;
Pw = 10.^(P/10)/1000;
np = numel(P);
T = c.metrics.turbulence;
nt = size(T,1);
[~,design] = min(abs(P-c.metrics.designPowerDBm));
[nz,nw] = size(C);
keys = {'direct','steering','flat_top','robust_aber','robust_op','max_power'};
labels = {'Direct (no RIS)','Steering-only','Flat top','Robust ABER','Robust OP','Max mean power'};
list = cell(nz,nw);
candidateRows = cell(nz,nw);
validationRows = cell(nz,nw);
rankRows = cell(nz,nw);
for iw = 1:nw
    for iz = 1:nz
        q = C{iz,iw};
        E = q.candidates;
        K = numel(E);
        family = string({E.family});
        spots = [E.spot];
        feasible = [spots.minorRadiusM]>=c.study.minimumWidthM & ...
            [spots.equivalentRadiusM]>=c.study.minimumWidthM & ...
            [spots.axisRatio]<=c.study.maximumAxisRatio & [spots.r50M]>=0.02;
        feasible(family=="steering") = true;
        psi = double(q.psi);
        accepted = psi<=c.study.fovHalfAngleRad;
        curves = cell(1,nt); choices = cell(1,nt); ranks = cell(1,nt);
        fprintf('[evaluate] U = %g m/s, z = %g m\n',q.wind,q.distanceM);
        for t = 1:nt
            [A,O,M] = metrics(c,E,'gainField',psi,Pw,T(t,:));
            d = rd_link_metrics(c,q.direct.gain,psi,Pw,T{t,1},T{t,2});
            ch = rd_select(A,O,M,family,feasible,c.metrics.policy,design);
            curves{t} = assemble(A,O,M,d,ch);
            choices{t} = int16(ch);
            ranks{t} = rd_rank(curves{t},P,c);
            if t==1
                [AM,OM,MM] = metrics(c,E,'gainMC',psi,Pw,T(t,:));
                mcCurves = assemble(AM,OM,MM,d,ch);
                main = struct('A',A(:,design),'O',O(:,design),'AM',AM(:,design),'OM',OM(:,design));
            end
        end
        ch = choices{1};
        picks = [0,ch(2,design,1),ch(3,design,1),ch(4,design,1),ch(4,design,2),ch(5,design,1)];
        maps = struct('key',keys,'label',labels,'candidate',num2cell(picks), ...
            'x',[],'y',[],'field',[],'mc',[],'spotCenter',[],'spotCovariance',[], ...
            'meanGainField',[],'meanGainMC',[]);
        h = q.direct.mapHalfWidthM;
        bins = size(q.direct.map,1);
        axisM = ((1:bins)-0.5)/bins*2*h-h;
        maps(1).x = q.targetXY(1)+axisM; maps(1).y = q.targetXY(2)+axisM;
        maps(1).field = q.direct.map; maps(1).mc = q.direct.map;
        maps(1).meanGainField = mean(double(q.direct.gain(accepted)));
        maps(1).meanGainMC = maps(1).meanGainField;
        for j = 2:numel(maps)
            e = E(picks(j));
            maps(j).x = e.x; maps(j).y = e.y; maps(j).field = e.mapField; maps(j).mc = e.mapMC;
            maps(j).spotCenter = e.spot.centerXY; maps(j).spotCovariance = e.spot.covarianceM2;
            maps(j).meanGainField = mean(double(e.gainField(accepted)));
            maps(j).meanGainMC = mean(double(e.gainMC(accepted)));
        end
        r = struct();
        r.distanceM = q.distanceM; r.wind = q.wind; r.Z = q.Z;
        r.targetXY = q.targetXY; r.rxRay = q.rxRay; r.q0 = q.q0; r.pose = q.pose;
        r.fovOutageExact = q.fovOutageExact;
        r.family = family; r.width = [E.width]; r.feasible = feasible;
        r.curves = curves; r.mcCurves = mcCurves; r.choice = choices; r.rank = ranks;
        r.maps = maps;
        r.direct = rmfield(q.direct,{'gain','map'});
        r.direct.meanGain = maps(1).meanGainField;
        list{iz,iw} = r;

        gF = arrayfun(@(e) mean(double(e.gainField(accepted))),E);
        gM = arrayfun(@(e) mean(double(e.gainMC(accepted))),E);
        candidateRows{iz,iw} = table(repmat(q.wind,K,1),repmat(q.distanceM,K,1),(1:K).', ...
            family.',[E.width].',feasible.',[spots.equivalentRadiusM].',[spots.minorRadiusM].', ...
            [spots.axisRatio].',[spots.r50M].',[spots.r90M].',[E.saturatedPowerFraction].', ...
            [E.maxCommandDeg].',gF.',gM.',(gM./gF-1).',[E.outsideMapFraction].', ...
            ([E.beamletTruncatedPowerW]./[E.pupilPowerW]).', ...
            arrayfun(@(e) e.coherenceRingPowerW(1)/max(e.mapPowerW,realmin),E).', ...
            [E.mapPowerW].',[E.mcMapPowerW].',main.A,main.AM,main.O,main.OM,[E.elapsedSeconds].', ...
            'VariableNames',{'windMps','distanceM','design','family','targetWidthM','feasible', ...
            'equivalentRadiusM','minorRadiusM','axisRatio','r50M','r90M','saturatedPowerFraction', ...
            'maxCommandDeg','meanGainField','meanGainMC','relativeDifferenceMC','outsideMapFraction', ...
            'beamletTruncatedFraction','ring1CrossPowerFraction','mapPowerFieldW','mapPowerMCW', ...
            'aberFieldAtDesign','aberMCAtDesign','opFieldAtDesign','opMCAtDesign','elapsedSeconds'});
        ids = picks(2:end).';
        validationRows{iz,iw} = table(repmat(q.wind,5,1),repmat(q.distanceM,5,1), ...
            string(keys(2:end)).',ids,gF(ids).',gM(ids).',(gM(ids)./gF(ids)-1).', ...
            [main.A(ids(1:3));main.O(ids(4));main.A(ids(5))], ...
            [main.AM(ids(1:3));main.OM(ids(4));main.AM(ids(5))], ...
            'VariableNames',{'windMps','distanceM','scheme','design','meanGainField','meanGainMC', ...
            'relativeDifferenceMC','riskFieldAtDesign','riskMCAtDesign'});
        rows = cell(nt,1);
        for t = 1:nt
            parts = cell(2,1);
            for k = 1:2
                name = {'aber','op'};
                rk = ranks{t}.(name{k});
                parts{k} = table(repmat(q.wind,5,1),repmat(q.distanceM,5,1),repmat(t,5,1), ...
                    repmat(string(T{t,1}),5,1),repmat(string(T{t,2}),5,1),repmat(string(upper(name{k})),5,1), ...
                    (1:5).',rk.requiredDBm(:),ismember((1:5).',rk.winners),repmat(rk.clear,5,1), ...
                    repmat(rk.gainOverSteeringDB,5,1),repmat(rk.gainOverDirectDB,5,1), ...
                    'VariableNames',{'windMps','distanceM','turbulenceId','ocean','atmosphere', ...
                    'metric','scheme','requiredDBm','winner','clear','bestGainOverSteeringDB', ...
                    'bestGainOverDirectDB'});
            end
            rows{t} = vertcat(parts{:});
        end
        rankRows{iz,iw} = vertcat(rows{:});
    end
end
R = struct();
R.version = c.version;
R.created = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
R.config = c;
R.powerDBm = P;
R.designPowerDBm = P(design);
R.distancesM = c.scene.distancesM;
R.winds = c.sea.windsMps;
R.mainIndex = [find(R.distancesM==c.main.distanceM,1),find(R.winds==c.main.windMps,1)];
R.turbulence = T;
R.schemeNames = {'Direct (no RIS)','Steering-only','Flat top','Robust ABER/OP','Max mean power'};
R.schemeShort = {'Dir','Str','Flat','Rob','MaxP'};
R.cond = reshape([list{:}],nz,nw);
R.candidateTable = vertcat(candidateRows{:});
R.validationTable = vertcat(validationRows{:});
R.rankTable = vertcat(rankRows{:});
schemes = string(R.schemeNames).';
R.rankTable.schemeName = schemes(R.rankTable.scheme);
end

function [A,O,M] = metrics(c,E,field,psi,Pw,row)
K = numel(E);
A = zeros(K,numel(Pw)); O = A; M = A;
for k = 1:K
    m = rd_link_metrics(c,E(k).(field),psi,Pw,row{1},row{2});
    A(k,:) = m.aber; O(k,:) = m.op; M(k,:) = m.meanPowerW;
end
end

function cv = assemble(A,O,M,d,ch)
np = size(A,2);
cv = struct('aber',zeros(5,np),'op',zeros(5,np),'meanPowerW',zeros(5,np));
cv.aber(1,:) = d.aber; cv.op(1,:) = d.op; cv.meanPowerW(1,:) = d.meanPowerW;
for s = 2:5
    ia = sub2ind(size(A),ch(s,:,1),1:np);
    io = sub2ind(size(A),ch(s,:,2),1:np);
    cv.aber(s,:) = A(ia); cv.op(s,:) = O(io); cv.meanPowerW(s,:) = M(ia);
end
end
