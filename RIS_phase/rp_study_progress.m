function rp_study_progress(c,stage,completed,total,status)
%RP_STUDY_PROGRESS Publish a small atomic progress record during long stages.
if nargin<3, completed=0; end
if nargin<4, total=0; end
if nargin<5, status='running'; end
value=struct('status',status,'stage',stage,'completed',completed,'total',total, ...
    'time',char(datetime('now','TimeZone','Asia/Shanghai','Format','yyyy-MM-dd HH:mm:ss')));
temporary=[tempname(c.output),'.json'];
fid=fopen(temporary,'w');
assert(fid>=0,'rp:Progress','Cannot write progress record.');
fprintf(fid,'%s',jsonencode(value)); fclose(fid);
[ok,message]=movefile(temporary,fullfile(c.output,'progress.json'),'f');
assert(ok,'rp:Progress','Cannot publish progress: %s',message);
end
