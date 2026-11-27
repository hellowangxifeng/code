root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; q=s.cases{1}; g=q.geometry; a=q.illumination;
ids=a.activeIds; ref=a.referenceXY(ids,:);
incoming=rp_layer_path(ref-g.sourceLocal(1:2).',[-g.sourceLocal(3),g.hIn],[c.optics.nWater,c.optics.nGlass],c.numerics.bisectionIterations);
qin=incoming.q;
B=c.optics.nAir*sind(c.device.maxAngleDeg);
fprintf('B=%.6f\n',B);
fprintf('qin min/max x: %.4f %.4f  y: %.4f %.4f\n',min(qin(:,1)),max(qin(:,1)),min(qin(:,2)),max(qin(:,2)));
lo=max(qin)-B; hi=min(qin)+B;
fprintf('feasible box x [%.4f %.4f], y [%.4f %.4f]\n',lo(1),hi(1),lo(2),hi(2));
% grid search q target in box, |q|<1; compute min slack (distance to box edges and 1-|q|)
qx=linspace(max(lo(1),-1),min(hi(1),1),401); qy=linspace(max(lo(2),-1),min(hi(2),1),401);
[QX,QY]=meshgrid(qx,qy); R=hypot(QX,QY); feas=R<1;
slackBox=min(min(QX-lo(1),hi(1)-QX),min(QY-lo(2),hi(2)-QY));
slackDisk=1-R;
slack=min(slackBox,slackDisk); slack(~feas)=-Inf;
[mx,ix]=max(slack(:)); fprintf('max slack %.4f at q=[%.4f %.4f], angle %.3f deg\n',mx,QX(ix),QY(ix),asind(hypot(QX(ix),QY(ix))));
% minimal |q| for slack >= 0.02, 0.04, 0.06, 0.08
for sl=[0 0.02 0.04 0.06 0.08]
  mask=feas & slack>=sl;
  if ~any(mask(:)), fprintf('slack %.2f infeasible\n',sl); continue; end
  rr=R; rr(~mask)=Inf; [mn,j]=min(rr(:));
  fprintf('slack>=%.2f: min angle %.3f deg at q=[%.4f %.4f]\n',sl,asind(mn),QX(j),QY(j));
end
% saturated mean output
bcur=[-B B]; qout=qin+bcur; d=hypot(qout(:,1),qout(:,2));
fprintf('current qout angle: mean %.3f std %.3f min %.3f max %.3f deg; spread(95-5) %.3f\n',mean(asind(d)),std(asind(d)),min(asind(d)),max(asind(d)),diff(prctile(asind(d),[5 95])));
fprintf('current qout mean q [%.4f %.4f], norm %.4f, angle %.3f\n',mean(qout(:,1)),mean(qout(:,2)),norm(mean(qout)),asind(norm(mean(qout))));
