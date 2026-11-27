function folders = wrd_setup()
%WRD_SETUP Put RIS_dynamic and the read-only project modules on the path.
%   folders = wrd_setup() returns the resolved folder names. The modules in
%   W_RIS, RIS_phase, sea_wave and yimiao/fig3-7 are only CALLED; nothing in
%   those folders is modified by RIS_dynamic.
here = fileparts(mfilename('fullpath'));   % ...\code\RIS_dynamic\algorithm
dynamic = fileparts(here);                  % ...\code\RIS_dynamic
root = fileparts(dynamic);                  % ...\code
% The HG sampler / water table folder is called yimiao in this repository
% and 3yimiao on some older copies; accept either.
mcCandidates = {fullfile(root,'yimiao','fig3-7'),fullfile(root,'3yimiao','fig3-7')};
mcExisting = mcCandidates(cellfun(@isfolder,mcCandidates));
mc = mcCandidates{1};
if ~isempty(mcExisting), mc = mcExisting{1}; end
folders = struct('algorithm',here,'plots',fullfile(dynamic,'plots'), ...
    'tests',fullfile(dynamic,'tests'),'dynamic',dynamic,'root',root, ...
    'wris',fullfile(root,'W_RIS'),'risPhase',fullfile(root,'RIS_phase'), ...
    'seaWave',fullfile(root,'sea_wave'),'mc',mc,'mcCandidates',{mcExisting});
required = {folders.wris,folders.risPhase,folders.seaWave,folders.mc};
missing = required(~cellfun(@isfolder,required));
if ~isempty(missing)
    error('wrd:Setup',['RIS_dynamic must sit next to W_RIS, RIS_phase, sea_wave and yimiao ' ...
        '(or 3yimiao) inside the same "code" folder. Missing:\n  %s'],strjoin(missing,newline));
end
% addpath puts the FIRST argument first on the path: RIS_dynamic wins any
% name clash, followed by the original project modules.
addpath(here,folders.plots,folders.wris,folders.risPhase,folders.seaWave,folders.mc);
end
