function [g,mc,a]=rp_study_scene(c,z)
%RP_STUDY_SCENE Static wetted board; inherited TX angle, water and beam.
g=wris_reference(c);
g.sourceLocal=g.RRis.'*(g.pTx-g.pRis);
g.directionLocal=g.RRis.'*g.dTx;
g.hIn=c.ris.thickness*c.optics.screenFraction;
g.hOut=c.ris.thickness-g.hIn; g.zTop=c.ris.thickness;
g.rxPlaneZ=g.pRis(3)+c.ris.thickness+z;
g.rxDesired=[g.pRis(1:2);g.rxPlaneZ];
g.time=0; g.incidenceDeg=acosd(g.directionLocal(3));
mc=rp_mc(c,g); a=rp_illuminate(c,g,mc);
end
