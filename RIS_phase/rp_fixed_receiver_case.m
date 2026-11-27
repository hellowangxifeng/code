function b=rp_fixed_receiver_case(c,d,family,width,margin)
%RP_FIXED_RECEIVER_CASE Full loaded phase propagation onto the common RX region.
rp_setup(); folder=fullfile(c.output,'fixed_receiver');
if ~isfolder(folder), mkdir(folder); end
p=rp_fixed_receiver_phase(c,d,family,width,margin);
target=p.targetCenterXY;
file=fullfile(folder,sprintf('%s_w%03d.mat',family,round(1000*width)));
digest=java.security.MessageDigest.getInstance('SHA-256');
key=[rp_stage2_optical_signature(c),jsonencode(p),jsonencode(c.study), ...
    fileread([mfilename('fullpath') '.m']),fileread(which('rp_fixed_receiver_phase')), ...
    fileread(which('rp_aperture_channel')),fileread(which('rp_link_metrics'))];
digest.update(unicode2native(key,'UTF-8'));
signature=lower(reshape(dec2hex(typecast(digest.digest(),'uint8'),2).',1,[]));
if isfile(file)
    old=load(file,'b','signature');
    if strcmp(old.signature,signature), b=old.b; return; end
end
r=rp_trace(c,d.geometry,d.illumination,p);
native=rp_receiver_window(c,d.geometry,d.illumination,p,r);
% Cover the field AND >5-sigma position errors around the one fixed receiver.
p.centerXY=target; p.receiverHalfWidthM=max(.50,abs(native.centerXY-target)+native.halfWidthM);
w=rp_wave(c,d.geometry,d.illumination,d.mc,p,r); spot=rp_spot_stats(w);
h=rp_aperture_channel(c,d.geometry,w,target,rp_link_samples(c));
m=rp_link_metrics(c,h,c.study.designPowerW,64);
b=struct('family',family,'targetWidthM',width,'margin',margin,'phase',p, ...
    'wave',w,'spot',spot,'channel',h,'designMetrics',m,'receiverCenterXY',target, ...
    'receiverDirection',h.nominalDirection);
save(file,'b','signature','-v7.3');
fprintf('Fixed receiver %s %.1f cm: ABER %.6g OP %.6g power %.6g mW\n', ...
    family,100*width,m.aber,m.op,1000*m.meanReceivedPowerW);
end
