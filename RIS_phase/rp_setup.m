function rp_setup()
%RP_SETUP Add the existing motion models and this simulation to the path.
root = fileparts(mfilename('fullpath'));
parent = fileparts(root);
folders={root,fullfile(root,'plots'),fullfile(parent,'W_RIS'),fullfile(parent,'sea_wave'), ...
    fullfile(parent,'yimiao','fig3-7'),fullfile(parent,'3yimiao','fig3-7')};
addpath(folders{cellfun(@isfolder,folders)});
end
