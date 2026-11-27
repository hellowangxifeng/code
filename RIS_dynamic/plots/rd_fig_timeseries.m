function files = rd_fig_timeseries(R,o)
%RD_FIG_TIMESERIES Fig4: one rough-sea record frame by frame.
%   (a) jitter-averaged gain of every link (no fading), dB;
%   (b) frame OP at o.fixedPowerDBm (jitter, FOV and EGG x GG inside a frame);
%   (c) board tilt, sea slope at the direct-beam crossing, and the fraction
%       of ballistic direct-link power lost to total internal reflection.
files = strings(0,1);
windy = find(R.winds>0);
if isempty(windy), return; end
iw = windy(rd_fig_index(R.winds(windy),o.timeWindMps));
iz = rd_fig_index(R.distancesM,o.timeDistanceM);
F = R.main.frames{iw};
t = F.t;
style = rd_fig_style();
keys = R.schemes.keys;
P = R.powerDBm;
[fig,layout] = rd_fig_new(3,1,o,5.6,1.75);
handles = gobjects(1,numel(keys));
tLimits = [t(1) max(t(end),t(1)+eps)];

ax = nexttile(layout); hold(ax,'on');
for k = 1:numel(keys)
    y = F.meanGainFrames(:,iz,k);
    if ~any(isfinite(y)), continue; end
    handles(k) = plot(ax,t,10*log10(max(y,1e-15)),'Color',style.colors(k,:), ...
        'LineStyle',style.lines{k},'LineWidth',style.widths(k));
end
ylabel(ax,'Gain (dB)'); grid(ax,'on'); xlim(ax,tLimits);
text(ax,0.01,0.95,sprintf('(a) z = %g m, U = %g m/s',R.distancesM(iz),R.winds(iw)), ...
    'Units','normalized','VerticalAlignment','top','FontName',o.fontName, ...
    'FontSize',o.fontSize,'BackgroundColor','w','Margin',1);

ax = nexttile(layout); hold(ax,'on');
for k = 1:numel(keys)
    frames = double(F.opFrames(:,:,iz,k));
    if ~any(isfinite(frames(:))), continue; end
    y = arrayfun(@(f) rd_fig_at(frames(f,:),P,o.fixedPowerDBm),1:size(frames,1));
    y(y<=0) = NaN;
    plot(ax,t,y,'Color',style.colors(k,:),'LineStyle',style.lines{k},'LineWidth',style.widths(k));
end
set(ax,'YScale','log','YLim',o.yLimits); grid(ax,'on'); xlim(ax,tLimits);
ylabel(ax,'Frame OP');
text(ax,0.01,0.95,sprintf('(b) P_t = %g dBm',o.fixedPowerDBm),'Units','normalized', ...
    'VerticalAlignment','top','FontName',o.fontName,'FontSize',o.fontSize, ...
    'BackgroundColor','w','Margin',1);

ax = nexttile(layout); hold(ax,'on');
yyaxis(ax,'left');
h1 = plot(ax,t,F.tiltDeg,'k-','LineWidth',1.1);
h2 = plot(ax,t,F.slopeDeg,'k--','LineWidth',1.0);
ylabel(ax,'Angle (deg)');
ax.YColor = [0 0 0];
yyaxis(ax,'right');
h3 = plot(ax,t,100*F.ballisticTirFraction,'-','Color',[.83 .16 .16],'LineWidth',1.1);
ylabel(ax,'Direct TIR (%)');
ax.YColor = [.83 .16 .16];
xlim(ax,tLimits); grid(ax,'on');
xlabel(ax,'t (s)');
legend(ax,[h1 h2 h3],{'RIS tilt','sea slope at beam','direct ballistic TIR'}, ...
    'Location','northeast','Box','off','FontName',o.fontName,'FontSize',o.fontSize-1, ...
    'Orientation','horizontal');
text(ax,0.01,0.95,'(c)','Units','normalized','VerticalAlignment','top', ...
    'FontName',o.fontName,'FontSize',o.fontSize);

rd_fig_legend(layout,handles,R.schemes.names,o);
files = rd_fig_save(fig,'fig4_time_series',o);
end
