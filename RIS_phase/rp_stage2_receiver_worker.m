function rp_stage2_receiver_worker(s)
%RP_STAGE2_RECEIVER_WORKER Timer-safe wrapper with durable failure evidence.
try
    errorFile=fullfile(s.config.output,'receiver_error.txt');
    if isfile(errorFile), delete(errorFile); end
    rp_stage2_receiver(s);
catch exception
    fid=fopen(fullfile(s.config.output,'receiver_error.txt'),'w','n','UTF-8');
    if fid>=0, fprintf(fid,'%s',getReport(exception,'extended','hyperlinks','off')); fclose(fid); end
    rp_study_progress(s.config,'receiver_scans_failed',0,3,'failed');
end
end
