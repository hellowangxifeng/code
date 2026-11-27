function rd_tables(R,folder)
%RD_TABLES CSV exports of the light result (one row per design / scheme / point).
%   designs.csv     every RIS design of every condition: spot, LC program,
%                   field and MC mean gains, numerical diagnostics
%   validation.csv  the designs picked at the design power: field vs MC
%   ranking.csv     required P_t to reach the ABER / OP targets, winners
%   curves.csv      the five scheme curves of the main fading row
writetable(R.candidateTable,fullfile(folder,'designs.csv'));
writetable(R.validationTable,fullfile(folder,'validation.csv'));
writetable(R.rankTable,fullfile(folder,'ranking.csv'));
[nz,nw] = size(R.cond);
P = R.powerDBm(:);
np = numel(P);
parts = cell(nz,nw);
for iw = 1:nw
    for iz = 1:nz
        q = R.cond(iz,iw);
        cv = q.curves{1};
        mc = q.mcCurves;
        ch = double(q.choice{1});
        block = cell(5,1);
        for s = 1:5
            block{s} = table(repmat(q.wind,np,1),repmat(q.distanceM,np,1),repmat(s,np,1), ...
                repmat(string(R.schemeNames{s}),np,1),P,cv.aber(s,:).',cv.op(s,:).', ...
                cv.meanPowerW(s,:).',mc.aber(s,:).',mc.op(s,:).',ch(s,:,1).',ch(s,:,2).', ...
                'VariableNames',{'windMps','distanceM','scheme','schemeName','transmitDBm', ...
                'aber','op','meanReceivedPowerW','aberMC','opMC','designForABER','designForOP'});
        end
        parts{iz,iw} = vertcat(block{:});
    end
end
writetable(vertcat(parts{:}),fullfile(folder,'curves.csv'));
end
