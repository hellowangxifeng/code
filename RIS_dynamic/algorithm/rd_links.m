function S = rd_links(c,designs,seaAware,samples)
%RD_LINKS Assemble the seven compared links and their fixed receivers.
%   S = rd_links(c,designs,seaAware,samples): designs is one RD_DESIGN result
%   per receiver height. seaAware(iz) is the bank index chosen by RD_TRAIN
%   for this wind ([] for the calm sea -> Robust ABER is used). Designs made
%   from explicit specs (incidence sweep) already carry a 'sea_aware' entry.
sch = rd_schemes();
K = numel(sch.keys);
nz = numel(designs);
S = struct('z',cellfun(@(d) d.distanceM,designs),'Z',cellfun(@(d) d.Z,designs), ...
    'g0',designs{1}.g,'crossXY',designs{1}.crossXY,'links',{cell(nz,K)},'signature',[]);
signature = repmat(struct('key','','family','','width',NaN,'margin',NaN, ...
    'center',[NaN NaN]),nz,K);
for iz = 1:nz
    d = designs{iz};
    for k = 1:K
        key = sch.keys{k};
        if strcmp(key,'direct')
            e = struct('family','direct','width',0,'margin',0,'p',[], ...
                'center',d.direct.center,'valid',true);
            origin = d.originDirect;
            halfWidth = c.maps.halfWidthDirectM;
        else
            e = pick(d,key,seaAware,iz);
            origin = d.originRis;
            halfWidth = c.maps.halfWidthRisM;
        end
        available = ~isempty(e) && e.valid;
        center = [NaN NaN];
        p = [];
        if available
            center = e.center;
            p = e.p;
        end
        S.links{iz,k} = struct('key',key,'available',available,'p',p, ...
            'rx',rd_receiver(c,origin,center,d.Z,samples),'mapHalfWidthM',halfWidth);
        signature(iz,k).key = key;
        if ~isempty(e)
            signature(iz,k).family = e.family;
            signature(iz,k).width = e.width;
            signature(iz,k).margin = e.margin;
        end
        signature(iz,k).center = center;
    end
end
S.signature = signature;
end

function e = pick(d,key,seaAware,iz)
id = find(strcmp({d.ris.key},key),1);
if ~isempty(id)
    e = d.ris(id);
    return
end
e = [];
if ~strcmp(key,'sea_aware'), return; end
if ~isempty(seaAware) && iz<=numel(seaAware) && isfinite(seaAware(iz))
    e = d.bank(seaAware(iz));
else
    id = find(strcmp({d.ris.key},'robust_aber'),1);
    if ~isempty(id), e = d.ris(id); end
end
end
