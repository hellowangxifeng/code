root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase'),fullfile(root,'FOWT_RIS_eval','src')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; qc=s.cases{1}; I=load(fullfile(root,'FOWT_RIS_eval','results','improved_nominal.mat'));
for which=1:2
  if which==1, g=qc.geometry; a=qc.illumination; p=qc.phase; tag='current'; else, g=qc.geometry; a=I.p; p=I.p; tag='improved'; end
  % build reference-ray landing
  ref=p.referenceXY; q=p.qInAtReference+p.gradientOverK0; nA=c.optics.nAir; nG=c.optics.nGlass;
  gk=sqrt(nG^2-sum(q.^2,2)); topB=ref+c.ris.thickness*0.5*q./gk;
  topW=[topB,repmat(g.zTop,size(q,1),1)]*g.RRis.'+g.pRis.';
  d=[q,sqrt(nA^2-sum(q.^2,2))]/nA*g.RRis.';
  dist=(g.rxPlaneZ-topW(:,3))./d(:,3); land=topW+dist.*d;
  % weights: unit direct power
  if which==1, ids=qc.illumination.activeIds; else, ids=p.unitIds; end
  P=qc.illumination.unitDirectPowerW(ids); P=P/sum(P);
  ctr=sum(land(:,1:2).*P,1);
  dev=land(:,1:2)-ctr;
  sd=sqrt(sum(P.*sum(dev.^2,2)));
  fprintf('%s reference-ray landing: centroid [%.4f %.4f], weighted rms radius %.5f m, max %.5f, aperture(5cm) fraction %.4f\n',tag,ctr,sd,max(vecnorm(dev,2,2)),sum(P(vecnorm(dev,2,2)<=0.05)));
  % direction spread weighted
  ang=atan2d(vecnorm(d(:,1:2),2,2),d(:,3)); mu=sum(P.*ang);
  fprintf('   weighted dir mean %.4f std %.4f deg; unweighted std %.4f\n',mu,sqrt(sum(P.*(ang-mu).^2)),std(ang));
end
