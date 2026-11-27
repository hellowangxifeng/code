function value = rd_cached(file,fingerprint,compute,resume)
%RD_CACHED Reuse a checkpoint only if its fingerprint (inputs + code) matches.
%   value = rd_cached(file,fingerprint,@() work(...),resume)
if resume && isfile(file)
    try
        saved = load(file,'value','fingerprint');
        if isfield(saved,'fingerprint') && isequaln(saved.fingerprint,fingerprint)
            [~,name] = fileparts(file);
            fprintf('  reusing checkpoint %s\n',name);
            value = saved.value;
            return
        end
        fprintf('  checkpoint %s is stale (inputs or code changed), recomputing\n',file);
    catch err
        fprintf('  checkpoint %s unreadable (%s), recomputing\n',file,err.message);
    end
end
value = compute();
folder = fileparts(file);
if ~isfolder(folder), mkdir(folder); end
save(file,'value','fingerprint','-v7.3');
end
