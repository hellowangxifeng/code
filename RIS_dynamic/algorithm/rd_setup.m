function folders = rd_setup()
%RD_SETUP Put RIS_dynamic and the read-only project modules on the path.
%   folders = rd_setup() returns the resolved folders. RIS_phase, W_RIS,
%   sea_wave and yimiao/fig3-7 are only CALLED; RIS_dynamic never writes
%   into them (all caches and results go to RIS_dynamic/results).
here = fileparts(mfilename('fullpath'));      % ...\code\RIS_dynamic\algorithm
dynamic = fileparts(here);                     % ...\code\RIS_dynamic
root = fileparts(dynamic);                     % ...\code
mcCandidates = {fullfile(root,'yimiao','fig3-7'),fullfile(root,'3yimiao','fig3-7')};
mcFound = mcCandidates(cellfun(@isfolder,mcCandidates));
folders = struct('algorithm',here,'plots',fullfile(dynamic,'plots'), ...
    'tests',fullfile(dynamic,'tests'),'dynamic',dynamic,'root',root, ...
    'risPhase',fullfile(root,'RIS_phase'),'wris',fullfile(root,'W_RIS'), ...
    'seaWave',fullfile(root,'sea_wave'),'mc','');
if ~isempty(mcFound), folders.mc = mcFound{1}; end
required = {folders.risPhase,folders.wris,folders.seaWave};
missing = required(~cellfun(@isfolder,required));
if isempty(folders.mc), missing{end+1} = mcCandidates{1}; end
if ~isempty(missing)
    error('rd:Setup',['RIS_dynamic must sit next to RIS_phase, W_RIS, sea_wave and ' ...
        'yimiao (or 3yimiao) in the same "code" folder. Missing:\n  %s'], ...
        strjoin(missing,newline));
end
addpath(folders.risPhase,fullfile(folders.risPhase,'plots'),folders.wris, ...
    folders.seaWave,folders.mc);
% addpath prepends: RIS_dynamic is added last so that it is searched first.
addpath(here,folders.plots);
end
