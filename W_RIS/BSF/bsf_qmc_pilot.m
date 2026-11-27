function audit = bsf_qmc_pilot(row,col,powerOfTwo)
%BSF_QMC_PILOT Measure independent-scramble variance before selecting production N.
arguments
    row (1,1) {mustBeInteger,mustBePositive,mustBeReal} = 3
    col (1,1) {mustBeInteger,mustBePositive,mustBeReal} = 3
    powerOfTwo (1,1) {mustBeInteger,mustBePositive,mustBeReal} = 18
end
root = bsf_setup();
data = load(fullfile(fileparts(root),'results','wris_results.mat'),'s');
c = data.s.config;
o = data.s.panels{row,col}.out;
g = o.geometry;
Ib = bsf_cell_average(@(x,y) bsf_ballistic(c,g,x,y),o.x,o.y,4);
maps = zeros([size(Ib),8]);
powers = zeros(8,1);
tic;
for j = 1:8
    trial = bsf_scatter_transport(c,g,Photons=2^powerOfTwo, ...
        Seed=82000+100*row+10*col+j,Sampling="sobol");
    maps(:,:,j) = trial.Is;
    powers(j) = trial.power;
end
se = std(maps,0,3)/sqrt(size(maps,3));
meanMap = mean(maps,3);
audit = struct('row',row,'col',col,'maps',maps,'mean',meanMap,'standardError',se, ...
    'powerMean',mean(powers),'powerSE',std(powers,0)/sqrt(numel(powers)), ...
    'relativeNoise',norm(se,'fro')/norm(Ib+meanMap,'fro'),'elapsedSeconds',toc);
audit.metrics = bsf_compare(Ib+meanMap,o.noTurbulence,o.x,o.y,o.cellArea,se,o.standardError/o.gain);
save(fullfile(root,sprintf('qmc_pilot_r%d_c%d.mat',row,col)),'audit');
fprintf('QMC r%d c%d, %d x 2^%d, noise %.5f%%, seconds %.2f\n', ...
    row,col,numel(powers),powerOfTwo,100*audit.relativeNoise,audit.elapsedSeconds);
disp(audit.metrics);
end
