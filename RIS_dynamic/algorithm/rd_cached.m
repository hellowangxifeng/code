function value = rd_cached(file,fingerprint,compute,resume)
%RD_CACHED Reuse a checkpoint only if its fingerprint (inputs + code) matches.
%   value = rd_cached(file,fingerprint,@() work(...),resume)
if rd_cache_valid(file,fingerprint,resume)
    saved = load(file,'value');
    [~,name] = fileparts(file);
    fprintf('  reusing checkpoint %s\n',name);
    value = saved.value;
    return
end
value = compute();
folder = fileparts(file);
if ~isfolder(folder), mkdir(folder); end
save(file,'value','fingerprint','-v7.3');
end
