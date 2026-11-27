function files = plot_wris_dynamic(resultFile)
%PLOT_WRIS_DYNAMIC Render saved results without running optical transport.
if nargin==0
    resultFile = fullfile(fileparts(mfilename('fullpath')),'results','wris_dynamic.mat');
end
loaded = load(resultFile,'results');
r = loaded.results;
folder = fileparts(resultFile);
if isempty(folder), folder = pwd; end
files = [plotPower(r,folder); plotMechanism(r,folder); plotMaps(r,folder)];
end

function paths = plotPower(results,folder)
shape = size(results.cases);
fig = newFigure([max(560,420*shape(2)),290*shape(1)+65]);
layout = tiledlayout(fig,shape(1),shape(2),'TileSpacing','compact','Padding','compact');
red = [0.74 0.29 0.19];
blue = [0.10 0.36 0.55];
tbl = results.summary;
positive = [tbl.directRxW;tbl.risRxW];
positive = positive(positive>0);
floorW = 1e-9;
if ~isempty(positive), floorW = min(positive)/3; end
topW = 1.25*max([floorW*50;tbl.directRxW+1.96*tbl.directRxSEW; ...
    tbl.risRxW+1.96*tbl.risRxSEW]);
for row = 1:shape(1)
    for col = 1:shape(2)
        item = results.cases{row,col};
        ax = nexttile(layout);
        hold(ax,'on');
        frames = item.frames;
        t = [frames.t];
        direct = [frames.direct];
        ris = [frames.ris];
        h1 = powerTrace(ax,t,direct,floorW,red);
        h2 = powerTrace(ax,t,ris,floorW,blue);
        h3 = plot(ax,NaN,NaN,'x','Color',[0.25 0.25 0.25],'MarkerSize',5);
        set(ax,'YScale','log');
        ylim(ax,1000*[floorW/1.6 topW]);
        xlabel(ax,'Elapsed time (s)');
        ylabel(ax,'RX aperture power (mW)');
        title(ax,caseTitle(item));
        grid(ax,'on');
        ax.YMinorGrid = 'off';
        if row==1 && col==1
            labels = {'Rough-sea direct','Moving RIS','No packet hit'};
            lgd = legend(ax,[h1,h2,h3],labels,'Orientation','horizontal');
            lgd.Layout.Tile = 'south';
            lgd.Box = 'off';
        end
    end
end
title(layout,'Same sea realization and fixed receiver');
paths = saveFigure(fig,folder,'paired_rx_power');
end

function h = powerTrace(ax,t,link,floorW,color)
value = [link.rxPowerW];
se = [link.powerSEW];
hasHit = [link.hitCount]>0;
shown = value;
shown(~hasHit) = NaN;
ci = 1.96*se;
h = errorbar(ax,t,shown*1e3,min(ci,max(0,value-floorW))*1e3,ci*1e3, ...
    '.-','Color',color,'LineWidth',1.15,'MarkerSize',9,'CapSize',0);
plot(ax,t(~hasHit),repmat(floorW*1e3,1,nnz(~hasHit)),'x', ...
    'Color',color,'MarkerSize',4,'HandleVisibility','off');
end

function paths = plotMechanism(results,folder)
shape = size(results.cases);
fig = newFigure([max(580,420*shape(2)),285*shape(1)+60]);
layout = tiledlayout(fig,shape(1),shape(2),'TileSpacing','compact','Padding','compact');
red = [0.74 0.29 0.19];
blue = [0.10 0.36 0.55];
for row = 1:shape(1)
    for col = 1:shape(2)
        item = results.cases{row,col};
        ax = nexttile(layout);
        f = item.frames;
        t = [f.t];
        direct = [f.direct];
        yyaxis(ax,'left');
        h1 = plot(ax,t,100*[direct.tirPowerW]./max(realmin,[direct.incidentPowerW]),'Color',red, ...
            'LineWidth',1.4);
        ylabel(ax,'Direct TIR power share (%)');
        ylim(ax,[0 100]);
        yyaxis(ax,'right');
        slopes = vertcat(f.slope);
        h2 = plot(ax,t,atand(slopes(:,1)),':','Color',blue,'LineWidth',1.5);
        ylabel(ax,'Sea slope along x (deg)');
        ax.YAxis(1).Color = red;
        ax.YAxis(2).Color = blue;
        xlabel(ax,'Elapsed time (s)');
        title(ax,caseTitle(item));
        grid(ax,'on');
        if row==1 && col==1
            lgd = legend(ax,[h1,h2],{'Power at sea surface undergoing TIR','Slope at nominal board centre'}, ...
                'Orientation','horizontal');
            lgd.Layout.Tile = 'south';
            lgd.Box = 'off';
        end
    end
end
title(layout,'Rough-sea refraction and total internal reflection');
paths = saveFigure(fig,folder,'rough_sea_tir');
end

function paths = plotMaps(results,folder)
% A preselected time and the strongest configured wind, shared by both links.
items = results.cases(:,end);
maps = cell(numel(items),2);
maxValue = 0;
times = zeros(numel(items),1);
for row = 1:numel(items)
    f = items{row}.frames;
    [~,k] = min(abs([f.t]-10));
    maps(row,:) = {f(k).direct,f(k).ris};
    maxValue = max([maxValue,max(f(k).direct.mapWm2,[],'all'),max(f(k).ris.mapWm2,[],'all')]);
    times(row) = f(k).t;
end
maxValue = max(maxValue,1e-7);
minimum = maxValue*1e-5;
fig = newFigure([870,350*numel(items)+45]);
layout = tiledlayout(fig,numel(items),2,'TileSpacing','compact','Padding','compact');
colormap(fig,parula(256));
linkNames = {'Rough-sea direct','Moving RIS'};
for row = 1:numel(items)
    for col = 1:2
        ax = nexttile(layout);
        link = maps{row,col};
        imagesc(ax,link.axisM,link.axisM,10*log10(max(minimum,link.mapWm2)));
        axis(ax,'image');
        ax.YDir = 'normal';
        clim(ax,10*log10([minimum maxValue]));
        hold(ax,'on');
        angle = linspace(0,2*pi,161);
        radius = results.config.receiver.radius;
        plot(ax,radius*cos(angle),radius*sin(angle),'w','LineWidth',1.2);
        xlabel(ax,'Detector plane u (m)');
        ylabel(ax,'Detector plane v (m)');
        title(ax,{caseTitle(items{row}),sprintf('%s, t = %.1f s',linkNames{col},times(row))});
    end
end
bar = colorbar(ax);
bar.Layout.Tile = 'east';
bar.Label.String = 'Irradiance before FOV gate (dB re 1 W m^{-2})';
title(layout,'Irradiance on the common detector plane');
paths = saveFigure(fig,folder,'paired_rx_maps');
end

function label = caseTitle(item)
label = sprintf('U = %.1f m/s, incidence = %.1f°',item.windMps, ...
    acosd(cosd(item.txRotationXYDeg(1))*cosd(item.txRotationXYDeg(2))));
end

function fig = newFigure(sizePx)
fig = figure('Visible','off','Color','w','Position',[100 100 sizePx]);
set(fig,'DefaultAxesFontName','Arial','DefaultAxesFontSize',10, ...
    'DefaultAxesLineWidth',0.75,'DefaultAxesBox','on');
end

function paths = saveFigure(fig,folder,name)
cleanup = onCleanup(@() close(fig));
paths = [string(fullfile(folder,[name '.png'])),string(fullfile(folder,[name '.fig']))];
exportgraphics(fig,paths(1),'Resolution',220);
savefig(fig,paths(2),'compact');
end
