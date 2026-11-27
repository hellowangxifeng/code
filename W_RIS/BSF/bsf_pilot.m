function pilot = bsf_pilot()
%BSF_PILOT Fast, reproducible numerical audit before the full 18-panel run.
root = bsf_setup();
data = load(fullfile(fileparts(root),'results','wris_results.mat'),'s');
c = data.s.config;
o = data.s.panels{1,1}.out;
g = o.geometry;
tic;
Ib = bsf_cell_average(@(x,y) bsf_ballistic(c,g,x,y),o.x,o.y,4);
scatter = bsf_scatter_transport(c,g,Photons=2e6,Seed=73101);
[X,Y] = ndgrid(o.x,o.y);
paper = bsf_paper(c,g,X,Y);
pilot = struct('elapsedSeconds',toc,'ballisticPower',sum(Ib,'all')*o.cellArea, ...
    'scatter',scatter,'paper',paper,'Ib',Ib);
pilot.transportMetrics = bsf_compare(Ib+scatter.Is,o.noTurbulence,o.x,o.y,o.cellArea, ...
    scatter.standardError,o.standardError/o.gain);
pilot.paperMetrics = bsf_compare(paper.IBSF,o.noTurbulence,o.x,o.y,o.cellArea, ...
    zeros(size(Ib)),o.standardError/o.gain);
save(fullfile(root,'pilot.mat'),'pilot');
disp(pilot.transportMetrics);
disp(pilot.paperMetrics);
fprintf('Pilot elapsed: %.2f s\n',pilot.elapsedSeconds);
end
