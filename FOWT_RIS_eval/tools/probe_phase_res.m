root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; q=s.cases{1}; g=q.geometry;
% choose the unit with max direct power
a=q.illumination; [~,im]=max(a.unitDirectPowerW); ref=a.referenceXY(im,:);
fprintf('unit ref [%.4f %.4f], max power %.4f\n',ref,a.unitDirectPowerW(im));
n=41; ex=linspace(-c.array.apertureM(1)/2,c.array.apertureM(1)/2,n);
[X,Y]=meshgrid(ex,ex); pts=[X(:),Y(:)]+ref;
f=rp_incident_field(c,g,pts);
ph=reshape(f.phaseRad,n,n);
% remove best plane fit (unwrap)
phu=unwrap(unwrap(ph,[],2),[],1);
A=[X(:) Y(:) ones(numel(X),1)]; coef=A\phu(:); res=phu(:)-A*coef;
fprintf('phase residual peak-to-peak %.2f rad = %.2f waves; rms %.2f rad\n',max(res)-min(res),(max(res)-min(res))/(2*pi),std(res));
% quadratic fit
A2=[X(:).^2 X(:).*Y(:) Y(:).^2 X(:) Y(:) ones(numel(X),1)]; coef2=A2\phu(:); res2=phu(:)-A2*coef2;
fprintf('quadratic residual p2p %.2f rad = %.2f waves; rms %.2f rad\n',max(res2)-min(res2),(max(res2)-min(res2))/(2*pi),std(res2));
