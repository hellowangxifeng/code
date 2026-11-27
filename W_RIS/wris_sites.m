function sites = wris_sites(c,wave)
%WRIS_SITES Four spatial sea positions in ONE realization at ONE fixed time.
% Positions refer to the initial RIS centre; subsequent free drift is retained.
direction=[sind(c.ris.windAngleDeg);cosd(c.ris.windAngleDeg)];
[omega,~,group]=unique(wave.omega);
energy=accumarray(group,abs(wave.coefficient).^2);
[~,peak]=max(energy);
wavelength=2*pi*c.ris.g/omega(peak)^2;
distance=linspace(-wavelength,wavelength,c.sea.siteSamples).';
points=c.ris.center.'+distance*direction.';
field=struct('wave',wave,'parameters',c.ris);
[z,gradient]=RIS_wave_field(field,points,c.sampling.comparisonTime);
slope=gradient*direction;
crest=find(z(2:end-1)>=z(1:end-2) & z(2:end-1)>z(3:end))+1;
trough=find(z(2:end-1)<=z(1:end-2) & z(2:end-1)<z(3:end))+1;
rising=find(z(1:end-1)<=0 & z(2:end)>0);
falling=find(z(1:end-1)>=0 & z(2:end)<0);
assert(~isempty(crest) && ~isempty(trough) && ~isempty(rising) && ~isempty(falling), ...
    'wris:Sites','Spatial transect lacks four sea positions; increase the transect span.');
[~,a]=max(z(crest)); [~,b]=min(z(trough));
[~,d]=max(slope(rising)); [~,e]=min(slope(falling));
ids=[crest(a),trough(b),rising(d),falling(e)];
labels={'Crest site','Trough site','Rising-slope site','Falling-slope site'};
sites=struct('label',labels,'center',cell(1,4),'elevation',cell(1,4),'slope',cell(1,4));
for j=1:4
    k=ids(j); coordinate=distance(k);
    if j>=3
        coordinate=fzero(@(q) RIS_wave_field(field,(c.ris.center+q*direction).', ...
            c.sampling.comparisonTime),distance(k:k+1));
    end
    sites(j).center=c.ris.center+coordinate*direction;
    [sites(j).elevation,g]=RIS_wave_field(field,sites(j).center.',c.sampling.comparisonTime);
    sites(j).slope=g*direction;
end
end
