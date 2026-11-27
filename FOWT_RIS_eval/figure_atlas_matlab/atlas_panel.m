function atlas_panel(ax, label, dy)
%ATLAS_PANEL Place (a)/(b) below the axes, reference-style.
if nargin < 3, dy = -0.20; end
text(ax, 0.5, dy, label, 'Units', 'normalized', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 14, 'Clipping', 'off');
end
