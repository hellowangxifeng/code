%RUN_BSF All 18 panels; resume saved batches until 0.1% RMS numerical SE.
bsf_setup();
try
    bsf_run_all(TargetRelativeSE=0.001,BatchPhotons=8e6,MaxPhotons=1024e6);
    bsf_validate_outputs();
catch failure
    root = fileparts(mfilename('fullpath'));
    fid = fopen(fullfile(root,'results','error_report.txt'),'w');
    fprintf(fid,'%s',getReport(failure,'extended','hyperlinks','off'));
    fclose(fid);
    rethrow(failure)
end
