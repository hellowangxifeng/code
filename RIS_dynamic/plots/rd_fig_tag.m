function rd_fig_tag(ax,text1,o,corner)
%RD_FIG_TAG Small condition label inside a free axes corner.
if nargin<4, corner = 'northeast'; end
switch corner
    case 'northeast', p = [0.97 0.96]; ha = 'right'; va = 'top';
    case 'northwest', p = [0.03 0.96]; ha = 'left'; va = 'top';
    case 'southwest', p = [0.03 0.04]; ha = 'left'; va = 'bottom';
    otherwise,        p = [0.97 0.04]; ha = 'right'; va = 'bottom';
end
text(ax,p(1),p(2),text1,'Units','normalized','HorizontalAlignment',ha,'VerticalAlignment',va, ...
    'FontName',o.fontName,'FontSize',o.fontSize,'BackgroundColor','w','Margin',1);
end
