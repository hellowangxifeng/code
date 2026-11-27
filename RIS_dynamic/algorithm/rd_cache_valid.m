function ok = rd_cache_valid(file,fingerprint,resume)
%RD_CACHE_VALID True if file holds a checkpoint with exactly this fingerprint.
% Only the small fingerprint variable is read, so the test is cheap.
ok = false;
if ~resume || ~isfile(file), return; end
try
    saved = load(file,'fingerprint');
    ok = isfield(saved,'fingerprint') && isequaln(saved.fingerprint,fingerprint);
    if ~ok
        fprintf('  checkpoint %s is stale (inputs or code changed), recomputing\n',file);
    end
catch err
    fprintf('  checkpoint %s unreadable (%s), recomputing\n',file,err.message);
end
end
