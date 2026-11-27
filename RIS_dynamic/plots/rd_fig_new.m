function [fig,layout] = rd_fig_new(rows,cols,o,tileWidthIn,tileHeightIn)
%RD_FIG_NEW White figure sized from the tile grid, with a compact tiled layout.
if nargin<4, tileWidthIn = o.tileWidthIn; end
if nargin<5, tileHeightIn = o.tileHeightIn; end
fig = figure('Visible',o.visible,'Color','w','Units','inches', ...
    'Position',[0.5 0.5 cols*tileWidthIn+0.7 rows*tileHeightIn+0.9]);
set(fig,'DefaultAxesFontName',o.fontName,'DefaultAxesFontSize',o.fontSize, ...
    'DefaultTextFontName',o.fontName,'DefaultAxesLineWidth',0.75,'DefaultAxesBox','on', ...
    'DefaultAxesTickDir','in');
layout = tiledlayout(fig,rows,cols,'TileSpacing','compact','Padding','compact');
end
