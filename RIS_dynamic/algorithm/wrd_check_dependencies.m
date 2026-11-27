function report = wrd_check_dependencies()
%WRD_CHECK_DEPENDENCIES Fail early, with a readable message, before a long run.
%   Verifies the MATLAB release and that every reused project function
%   resolves to the expected project folder. Only base MATLAB is required.
folders = wrd_setup();
try
    tooOld = isMATLABReleaseOlderThan("R2022a");
catch
    tooOld = true; % isMATLABReleaseOlderThan itself appeared in R2020b
end
if tooOld
    error('wrd:Release','RIS_dynamic needs MATLAB R2022a or newer (clim, tiledlayout, exportgraphics).');
end
needed = {
    'wris_config',folders.wris; 'wris_device',folders.wris; 'wris_mc',folders.wris; ...
    'wris_tx_motion',folders.wris; 'wris_reference',folders.wris; 'wris_intersect',folders.wris; ...
    'wris_metrics',folders.wris; ...
    'rp_config',folders.risPhase; 'rp_setup',folders.risPhase; 'rp_geometry',folders.risPhase; ...
    'rp_trace',folders.risPhase; 'rp_phase_profile',folders.risPhase; 'rp_fresnel',folders.risPhase; ...
    'rp_command',folders.risPhase; 'rp_layer_path',folders.risPhase; ...
    'RIS_float',folders.seaWave; 'RIS_wave_components',folders.seaWave; ...
    'RIS_wave_field',folders.seaWave; 'RIS_rotation',folders.seaWave; ...
    'RIS_flotation',folders.seaWave; 'RIS_radiation',folders.seaWave; ...
    'RIS_sidewall_load',folders.seaWave; 'get_default_RIS_params',folders.seaWave; ...
    'get_default_sea_params',folders.seaWave; ...
    'Table_abc',folders.mc; 'Hen_Green',folders.mc};
resolved = strings(size(needed,1),1);
problems = strings(0,1);
for j = 1:size(needed,1)
    location = which(needed{j,1});
    resolved(j) = string(location);
    expected = string(needed{j,2});
    if strcmp(needed{j,2},folders.mc), expected = string(folders.mcCandidates); end
    if isempty(location)
        problems(end+1) = sprintf('%s: not found (expected in %s)',needed{j,1},needed{j,2}); %#ok<AGROW>
    elseif ~any(strcmpi(string(fileparts(location)),expected))
        problems(end+1) = sprintf('%s: resolves to %s, expected %s', ...
            needed{j,1},location,needed{j,2}); %#ok<AGROW>
    end
end
if ~isempty(problems)
    error('wrd:Dependencies','Dependency check failed:\n  %s',strjoin(problems,newline+"  "));
end
report = table(string(needed(:,1)),resolved,'VariableNames',{'function','file'});
end
