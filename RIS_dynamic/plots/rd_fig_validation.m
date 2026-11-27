function files = rd_fig_validation(R,o)
%RD_FIG_VALIDATION Fig6: complex-field formula versus full Monte Carlo.
%   (a) In-FOV mean aperture gain <G | psi <= Theta> of EVERY RIS design of
%       every condition: complex field (x) versus MC (y); coloured markers
%       are the designs picked at the design power. Dashed: +-10 %.
%   (b) Main condition, irradiance averaged over the strip |y - y_RX| <= a
%       through the receiver: complex field (lines) versus MC (markers).
style = rd_fig_style();
T = R.candidateTable;
[fig,layout] = rd_fig_new(1,2,o,3.0,2.7);
ax = nexttile(layout);
hold(ax,'on');
ok = T.meanGainField>0 & T.meanGainMC>0;
plot(ax,T.meanGainField(ok),T.meanGainMC(ok),'.','Color',[.6 .6 .6],'MarkerSize',6, ...
    'DisplayName','all designs');
V = R.validationTable;
schemes = ["steering","flat_top","robust_aber","robust_op","max_power"];
color = [2 3 4 4 5];
marker = {'s','^','o','v','d'};
names = {'Steering-only','Flat top','Robust ABER','Robust OP','Max mean power'};
for j = 1:numel(schemes)
    use = V.scheme==schemes(j) & V.meanGainField>0 & V.meanGainMC>0;
    plot(ax,V.meanGainField(use),V.meanGainMC(use),marker{j},'Color',style.colors(color(j),:), ...
        'MarkerSize',5,'LineWidth',0.9,'DisplayName',names{j});
end
values = [T.meanGainField(ok);T.meanGainMC(ok)];
if isempty(values), values = [1e-6;1]; end
span = [10^floor(log10(min(values))),10^ceil(log10(max(values)))];
plot(ax,span,span,'k-','LineWidth',0.8,'HandleVisibility','off');
plot(ax,span,1.1*span,'k--','LineWidth',0.6,'HandleVisibility','off');
plot(ax,span,span/1.1,'k--','LineWidth',0.6,'HandleVisibility','off');
set(ax,'XScale','log','YScale','log','XLim',span,'YLim',span,'FontName',o.fontName, ...
    'FontSize',o.fontSize,'Box','on');
grid(ax,'on');
xlabel(ax,'complex field  \langleG | \psi\leq\Theta\rangle');
ylabel(ax,'Monte Carlo  \langleG | \psi\leq\Theta\rangle');
legend(ax,'Location','northwest','Box','off','FontName',o.fontName,'FontSize',o.fontSize-1);

iz = rd_fig_index(R.distancesM,o.mainDistanceM);
iw = rd_fig_index(R.winds,o.mainWindMps);
q = R.cond(iz,iw);
a = R.config.rx.apertureRadius;
ax = nexttile(layout);
hold(ax,'on');
keys = {q.maps.key};
colorOf = struct('direct',1,'steering',2,'flat_top',3,'robust_aber',4,'robust_op',4,'max_power',5);
for j = 1:numel(o.profileSchemes)
    k = find(strcmp(keys,o.profileSchemes{j}),1);
    if isempty(k), continue; end
    m = q.maps(k);
    rows = abs(m.y-q.targetXY(2))<=a;
    if ~any(rows), continue; end
    field = mean(double(m.field(rows,:)),1);
    mc = mean(double(m.mc(rows,:)),1);
    s = colorOf.(m.key);
    x = 100*(m.x-q.targetXY(1));
    field(field<=0) = NaN;
    mc(mc<=0) = NaN;
    plot(ax,x,field,'-','Color',style.colors(s,:),'LineWidth',o.lineWidth, ...
        'DisplayName',[m.label ', field']);
    step = max(1,round(numel(x)/40));
    plot(ax,x(1:step:end),mc(1:step:end),style.markers{max(s,2)},'Color',style.colors(s,:), ...
        'MarkerSize',4,'LineWidth',0.8,'DisplayName',[m.label ', MC']);
end
set(ax,'YScale','log','FontName',o.fontName,'FontSize',o.fontSize,'Box','on', ...
    'XLim',100*o.mapHalfWidthM*[-1 1]);
grid(ax,'on');
xlabel(ax,'x - x_{RX} (cm)');
ylabel(ax,'irradiance (W m^{-2})');
legend(ax,'Location','southoutside','Box','off','NumColumns',2,'FontName',o.fontName, ...
    'FontSize',o.fontSize-1);
rd_fig_tag(ax,sprintf('z = %g m, U = %g m/s',q.distanceM,q.wind),o,'northeast');
files = rd_fig_save(fig,'fig6_field_vs_mc',o);
end
