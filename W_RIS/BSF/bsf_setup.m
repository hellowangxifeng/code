function root = bsf_setup()
%BSF_SETUP Locate original modules without changing their source files.
root = fileparts(mfilename('fullpath'));
project = fileparts(root);
workspace = fileparts(project);
addpath(root, project, fullfile(workspace, 'sea_wave'), ...
    fullfile(workspace, '3yimiao', 'fig3-7'));
end
