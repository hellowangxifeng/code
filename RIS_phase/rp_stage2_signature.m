function key=rp_stage2_signature(c)
%RP_STAGE2_SIGNATURE Hash complete algorithm sources and current configuration.
digest=java.security.MessageDigest.getInstance('SHA-256'); cfg=c;
if isfield(cfg.study,'signature'), cfg.study=rmfield(cfg.study,'signature'); end
digest.update(unicode2native(jsonencode(cfg),'UTF-8'));
folders={c.root,fullfile(fileparts(c.root),'W_RIS')};
for folder=folders
    files=dir(fullfile(folder{1},'*.m'));
    for j=1:numel(files)
        digest.update(unicode2native(files(j).name,'UTF-8'));
        digest.update(unicode2native(fileread(fullfile(files(j).folder,files(j).name)),'UTF-8'));
    end
end
for name={'RIS_flotation','RIS_rotation'}
    digest.update(unicode2native(fileread(which(name{1})),'UTF-8'));
end
key=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
end
