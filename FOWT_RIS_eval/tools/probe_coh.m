root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'RIS_phase')); rp_setup();
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s; I=load(fullfile(root,'FOWT_RIS_eval','results','improved_nominal.mat'));
for which=1:2
  if which==1, w=s.cases{1}.wave; tag='current'; else, w=I.wave; tag='improved'; end
  coh=w.ballisticPixelMeanWm2; incoh=w.ballisticIncoherentWm2;
  fprintf('%s: max coh %.4f, max incoh %.4f, ratio %.2f; mean coh %.4f mean incoh %.4f ratio %.2f; sum coh %.4f sum incoh %.4f\n', ...
    tag,max(coh(:)),max(incoh(:)),max(coh(:))/max(incoh(:)),mean(coh(:)),mean(incoh(:)),mean(coh(:))/mean(incoh(:)),sum(coh(:))*w.cellAreaM2,sum(incoh(:))*w.cellAreaM2);
end
