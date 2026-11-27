root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
S=load(fullfile(root,'RIS_phase','results','validated_optics','ris_spots.mat'));
s=S.s;
for k=1:numel(s.cases)
  q=s.cases{k}; w=q.wave; u=q.illumination;
  fprintf('case %d: boardMC %.4f active %.4f direct %.4f | pupil %.4f inputDirect %.4f interfaceLoss %.4f | map %.4f diag %.4f scatter %.4f ballisticMC %.4f omitted %.2e | unresolved %.2e\n', ...
    k,w.ballisticBoardMCPowerW,u.activePowerW,u.activeDirectPowerW,w.pupilPowerW,w.inputDirectPowerW,w.interfaceLossW,w.mapPowerW,w.diagonalMapPowerW,w.scatteredMapPowerW,w.ballisticBoardMCPowerW,w.omittedBallisticPowerW,w.unresolvedQuadraturePowerFraction);
end
