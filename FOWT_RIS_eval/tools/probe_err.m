root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase'),fullfile(root,'FOWT_RIS_eval','src')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; c=s.config; q=s.cases{1}; g=q.geometry; a=q.illumination;
I=load(fullfile(root,'FOWT_RIS_eval','results','improved_nominal.mat'));
p=I.p; P=a.unitDirectPowerW(p.unitIds); P=P/sum(P);
u=I.uWorld; target.mode='direction'; target.u=u;
% recompute internals
ids=a.activeIds; ref=a.referenceXY(ids,:);
incoming=rp_layer_path(ref-g.sourceLocal(1:2).',[-g.sourceLocal(3),g.hIn],[c.optics.nWater,c.optics.nGlass],c.numerics.bisectionIterations);
qIn=incoming.q;
d=repmat(u, numel(ids),1); qT=c.optics.nAir*(d*g.RRis); qT=qT(:,1:2);
err=vecnorm(p.gradientOverK0-(qT-qIn),2,2);
fprintf('weighted rms b error %.5f; max %.5f\n',sqrt(sum(P.*err.^2)),max(err));
fprintf('power fraction with |b err|>0.002: %.5f\n',sum(P(err>0.002)));
% weighted direction spread
dref=p.referenceDirectionWorld;
ang=p.referenceAngleToWorldZDeg;
fprintf('power-weighted mean angle %.4f, std %.4f, p5 %.4f p95 %.4f\n',sum(P.*ang),sqrt(sum(P.*(ang-sum(P.*ang)).^2)),prctile(ang,5),prctile(ang,95));
% bins vs power
for lim=[0 0.001 0.01 0.1 1]
  m=P>lim; fprintf('P>%.3f: n=%d, wstd %.4f deg, weighted angle %.4f\n',lim,sum(m),sqrt(sum(P(m).*(ang(m)-sum(P(m).*ang(m)))^2)/sum(P(m))),sum(P(m).*ang(m))/sum(P(m)));
end
% check q target feasibility exactly
B=c.optics.nAir*sind(c.device.maxAngleDeg);
bad=any(abs(qT-qIn)>B+1e-12,2);
fprintf('units beyond box: %d, power frac %.6f\n',sum(bad),sum(P(bad)));
fprintf('max |qT| %.5f\n',max(vecnorm(qT,2,2)));
