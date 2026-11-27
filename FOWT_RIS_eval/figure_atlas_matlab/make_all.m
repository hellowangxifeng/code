%MAKE_ALL Generate the MATLAB journal-style figure atlas.
clear; clc; close all;
atlas_style;
S = load(fullfile('data','atlas_data.mat'));
fn = fieldnames(S);
for i=1:numel(fn)
    if isnumeric(S.(fn{i})) && ~isa(S.(fn{i}),'double')
        S.(fn{i}) = double(S.(fn{i}));
    end
end
fig01_flow();
fig02_slope(S);
fig03_slope(S);
fig05_link(S);
fig06_spots(S);
fig07_motion(S);
fig08_distance(S);
fig09_turb(S);
fig10_dynamic(S);
fig11_conditions(S);
fig14_centroid(S);
fig16_incidence(S);
fig17_aperture(S);
fig18_op(S);
fig19_aber(S);
fig20_distance_turb(S);
fig21_delay(S);
fig22_update(S);
fprintf('All MATLAB atlas figures generated.\n');
