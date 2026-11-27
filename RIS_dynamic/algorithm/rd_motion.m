function m = rd_motion(c,wind,seed)
%RD_MOTION Directional sea and floating-RIS six-DOF motion for one wind (cached).
%   One set of wave components (sea_wave/RIS_wave_components, Neumann
%   spectrum with cos^2 spreading) drives BOTH the board motion
%   (sea_wave/RIS_float) and, through m.wave, the rough surface crossed by
%   the direct link (RD_SURFACE). seed fixes the random wave phases.
cfg = c.ris;
cfg.windSpeed = wind;
cfg.windAngleDeg = c.sea.windDirectionDeg;
cfg.seed = seed;
cfg.plot = false;
cfg.wave = [];
file = fullfile(c.cache.folder,sprintf('motion_U%g_dir%g_seed%d.mat',wind, ...
    cfg.windAngleDeg,seed));
fingerprint = struct('cfg',cfg,'quadrature',c.sea.quadrature, ...
    'neumannConstant',c.sea.neumannConstant,'code',rd_code_hash({'RIS_float', ...
    'RIS_wave_components','RIS_flotation','RIS_radiation','RIS_rotation', ...
    'RIS_sidewall_load','get_default_RIS_params','get_default_sea_params'}));
m = rd_cached(file,fingerprint,@() solve(cfg,c.sea),c.run.resume);
end

function m = solve(cfg,sea)
fprintf('  solving RIS six-DOF motion: U = %g m/s, direction %g deg, seed %d\n', ...
    cfg.windSpeed,cfg.windAngleDeg,cfg.seed);
defaults = get_default_sea_params();
wave = RIS_wave_components(cfg,sea.quadrature);
wave.coefficient = wave.coefficient*sqrt(sea.neumannConstant/defaults.neumannConstant);
cfg.wave = wave;
m = RIS_float(cfg);
end
