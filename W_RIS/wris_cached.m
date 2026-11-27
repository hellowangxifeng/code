function value = wris_cached(c,kind,key,producer)
%WRIS_CACHED Cache physical computations without invalidating them for plot edits.
if ~c.cache.enabled
    value=producer();
    return
end
parent=fileparts(c.root);
switch kind
    case 'motion'
        names={'RIS_float.m','RIS_wave_components.m','RIS_wave_field.m','RIS_rotation.m', ...
            'RIS_radiation.m','RIS_sidewall_load.m','RIS_flotation.m', ...
            'get_default_RIS_params.m','get_default_sea_params.m'};
        files=cellfun(@(n) fullfile(parent,'sea_wave',n),names,'UniformOutput',false);
    case 'optics'
        files={fullfile(c.root,'wris_mc.m'),fullfile(c.root,'wris_intersect.m'), ...
            fullfile(c.root,'wris_metrics.m'),fullfile(parent,'3yimiao','fig3-7','Hen_Green.m')};
    otherwise
        error('wris:CacheKind','Unknown physical cache kind: %s.',kind);
end
digest=java.security.MessageDigest.getInstance('SHA-256');
digest.update(unicode2native(jsonencode(serializable(key)),'UTF-8'));
for j=1:numel(files)
    digest.update(unicode2native(fileread(files{j}),'UTF-8'));
end
hash=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
if ~isfolder(c.cache.folder), mkdir(c.cache.folder); end
file=fullfile(c.cache.folder,[kind,'_',hash,'.mat']);
if isfile(file)
    data=load(file,'value'); value=data.value;
else
    value=producer();
    save(file,'value','key','-v7.3');
end
end

function value=serializable(value)
% Preserve complex wave coefficients in JSON keys.
if isstruct(value)
    names=fieldnames(value);
    for j=1:numel(names), value.(names{j})=serializable(value.(names{j})); end
elseif iscell(value)
    value=cellfun(@serializable,value,'UniformOutput',false);
elseif isnumeric(value) && ~isreal(value)
    value=struct('realPart',real(value),'imaginaryPart',imag(value));
end
end
