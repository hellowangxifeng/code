function bsf_worker()
%BSF_WORKER Run precision study as a timer callback in the connected MATLAB.
root = fileparts(mfilename('fullpath'));
try
    bsf_run_all(TargetRelativeSE=0.001,BatchPhotons=8e6,MaxPhotons=1024e6);
catch failure
    fid = fopen(fullfile(root,'results','error_report.txt'),'w');
    fprintf(fid,'%s',getReport(failure,'extended','hyperlinks','off'));
    fclose(fid);
    fid = fopen(fullfile(root,'results','progress.json'),'w');
    fprintf(fid,'%s',jsonencode(struct('status','failed','message',failure.message)));
    fclose(fid);
    warning('bsf:Worker','BSF precision calculation failed: %s',failure.message);
end
end
