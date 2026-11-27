function wris_setup()
%WRIS_SETUP Add only required existing project modules.
root = fileparts(mfilename('fullpath'));
addpath(root,fullfile(fileparts(root),'sea_wave'),fullfile(fileparts(root),'3yimiao','fig3-7'));
end
