function s = wris_simulate(c)
%WRIS_SIMULATE Compute snapshots only; no graphics, reports or figure layout.
t=c.sampling.comparisonTime;
validateattributes(c.comparison.rotationXYDeg,{'numeric'},{'2d','ncols',2,'nonempty','real','finite'});
validateattributes(c.comparison.conditions,{'numeric'},{'2d','ncols',2,'nonempty','real','finite'});
assert(t>=0 && t<=c.ris.duration,'wris:Time','Snapshot time must lie inside the simulated interval.');
angles=c.comparison.rotationXYDeg;
% Reject invalid directions before an expensive hydrodynamic solve.
for row=1:size(angles,1)
    cc=c; cc.geometry.rotationXYDeg=angles(row,:); wris_reference(cc);
end
s.config=c; s.panels=cell(size(angles,1),size(c.comparison.conditions,1)+1);
site=struct('label','Initial centre','center',c.ris.center);
for row=1:size(angles,1)
    cc=c; cc.geometry.rotationXYDeg=angles(row,:);
    g=wris_geometry(cc,struct('t',t),struct(),1,'static');
    s.panels{row,1}=panel(cc,g,0,0,[],[],site);
end
for col=1:size(c.comparison.conditions,1)
    cc=c; cc.ris.windSpeed=c.comparison.conditions(col,1);
    cc.ris.windAngleDeg=c.comparison.conditions(col,2);
    [wave,~]=RIS_wave_components(cc.ris,c.sea.quadrature);
    defaults=get_default_sea_params();
    wave.coefficient=wave.coefficient*sqrt(c.sea.neumannConstant/defaults.neumannConstant);
    site=selectSite(cc,wave); cc.ris.center=site.center; cc.ris.wave=wave;
    fprintf('RIS motion: U=%g m/s, wind=%g deg, site=%s\n',cc.ris.windSpeed,cc.ris.windAngleDeg,site.label);
    m=wris_cached(cc,'motion',cc.ris,@() RIS_float(cc.ris));
    [~,k]=min(abs(m.t-t));
    assert(abs(m.t(k)-t)<1e-10,'wris:Time','Snapshot time must coincide with a saved motion sample.');
    for row=1:size(angles,1)
        cc.geometry.rotationXYDeg=angles(row,:);
        fprintf('Optics: [ax ay]=[%g %g] deg | U=%g | t=%g s\n',angles(row,:),cc.ris.windSpeed,t);
        tx=wris_tx_motion(cc,m.t); g=wris_geometry(cc,m,tx,k,'joint');
        s.panels{row,col+1}=panel(cc,g,cc.ris.windSpeed,cc.ris.windAngleDeg,m,tx,site);
    end
end
end

function site=selectSite(c,wave)
name=validatestring(c.comparison.site,{'center','crest','trough','rising','falling'});
if strcmp(name,'center')
    site=struct('label','Initial centre','center',c.ris.center);
else
    names={'crest','trough','rising','falling'};
    sites=wris_sites(c,wave); site=sites(strcmp(name,names));
end
end

function p=panel(c,g,wind,direction,motion,tx,site)
gain=1;
if ~isempty(motion)
    gains=wris_egg(c.turbulence.kind,numel(motion.t),c.turbulence.seed);
    [~,k]=min(abs(motion.t-g.time)); gain=gains(k);
end
key=struct('source',c.source,'water',c.water,'mc',c.mc,'grid',c.grid, ...
    'size',[c.ris.Lx c.ris.Ly],'geometry',g,'gain',gain);
out=wris_cached(c,'optics',key,@() wris_mc(c,g,c.mc.seed,gain));
ref=wris_reference(c);
p=struct('initialRotationXYDeg',c.geometry.rotationXYDeg,'initialIncidenceDeg',ref.thetaDeg, ...
    'windMps',wind,'windDirectionDeg',direction,'phase',site.label,'site',site, ...
    'motion',motion,'tx',tx,'out',out);
end
