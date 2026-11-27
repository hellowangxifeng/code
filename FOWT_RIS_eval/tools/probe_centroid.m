root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s;
for k=1:numel(s.cases)
  q=s.cases{k}; r=q.rays; f=q.footprint;
  ok=r.valid & isfinite(r.landing(:,1));
  w=r.powerW(ok); L=r.landing(ok,1:2);
  c=sum(L.*w,1)/sum(w);
  d=r.directionWorld(ok,:); dw=sum(d.*w,1)/sum(w); dw=dw/norm(dw);
  ang=acosd(dw(3));
  fc=f.centerXY(:).';
  fprintf('case %d ballistic centroid [%.4f %.4f] vs spot center [%.4f %.4f] diff %.4f m; mean dir angle %.3f vs %.3f deg; valid rays %d\n', ...
     k,c(1),c(2),fc(1),fc(2),norm(c-fc),ang,mean(q.phase.referenceAngleToWorldZDeg),sum(ok));
end
