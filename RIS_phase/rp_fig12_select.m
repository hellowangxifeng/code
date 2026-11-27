function best=rp_fig12_select(T)
%RP_FIG12_SELECT Independent family minima with explicit deterministic ties.
arguments
    T table
end
sigmaIds=unique(T.sigmaId); rows=cell(numel(sigmaIds)*5,1); count=0;
for k=sigmaIds.'
    for objective=1:5
        pool=T(T.sigmaId==k,:);
        if objective<=2, pool=pool(pool.family=="robust",:); end
        if objective>=4, pool=pool(pool.family=="flat_top",:); end
        pool.negativePower=-pool.meanReceivedPowerW;
        if ismember(objective,[1 4])
            pool=sortrows(pool,{'aber','op','negativePower','candidate'});
        elseif ismember(objective,[2 5])
            pool=sortrows(pool,{'op','aber','negativePower','candidate'});
        else
            pool=sortrows(pool,{'negativePower','aber','op','candidate'});
        end
        count=count+1; q=pool(1,:); q.negativePower=[];
        q.objectiveId=objective; rows{count}=q;
    end
end
best=vertcat(rows{:});
end
