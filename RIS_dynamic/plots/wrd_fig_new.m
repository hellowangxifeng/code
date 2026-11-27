function fig = wrd_fig_new(sizePx,opts)
%WRD_FIG_NEW White figure with the shared font settings.
%   Legends and colorbars set their fonts explicitly where they are created.
fig = figure('Visible',opts.visible,'Color','w','Units','pixels','Position',[80 80 sizePx]);
set(fig,'DefaultAxesFontName',opts.fontName,'DefaultAxesFontSize',opts.fontSize, ...
    'DefaultAxesLineWidth',0.75,'DefaultAxesBox','on','DefaultTextFontName',opts.fontName);
end
