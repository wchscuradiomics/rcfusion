%%
clear; clc;
load our.mat;

dsname = 'tst';
% x1~x5: rad-score, and cliniomics.Properties.VariableNames(param.fs.new.newsubset(2:end)-1)
% {'Rad-score';'Age';'Diabetes';'CAD';'NIHSS score'}
varnames = {'Rad-score','Age','Diabetes','CAD','NIHSS score'};
f = @(D) predict(nom,D);

switch dsname
  case 'all'
    explainer = shapley(f,TRN,'QueryPoints',[TRN;TST;TSTSCNO4H]); 
    cexplainer = cshapley(f,TRN,'QueryPoints',[TRN;TST;TSTSCNO4H]);
  case 'tst'
    explainer = shapley(f,TRN,'QueryPoints',TST);
    cexplainer = cshapley(f,TRN,'QueryPoints',TST); 
  case 'scno4h'
    explainer = shapley(f,TRN,'QueryPoints',TSTSCNO4H); 
    cexplainer = cshapley(f,TRN,'QueryPoints',TSTSCNO4H);
end

lblfontsize = 15; lblfontname = 'Helvetica';

%% Plot Dependence for each variables
for i=2:length(varnames)
  fig = cplotdependence(explainer,1,i,varnames,lblfontsize,lblfontname);
  svgfilename = ['shapley_' dsname '_dependence_ ' varnames{1} '_' varnames{i} '.svg'];
  saveas(fig,svgfilename,'svg');
  close(fig);
  trimSVGFile(svgfilename);
end

%% Plot mean(abs(shap)) for this multi-query-point shapley object
bars = plot(cexplainer);
fig = bars(1).Parent.Parent; ax = fig.CurrentAxes;
set(ax,'TickLabelInterpreter','none','FontName','Helvetica Narrow','FontSize',lblfontsize); % global setting for fig.CurrentAxes

h = findobj(fig, 'Type', 'Legend');
h.Position(4) = 0.15;
h.Position(2) = h.Position(2) + 0.05;

yls = yticklabels;
for i=1:length(yls), yls{i} = varnames{str2double(replace(yls{i},'x',''))}; end
yticklabels(yls);

legend boxoff;

ylabel('');
title('');
ax.Title.FontName = lblfontname;
ax.Title.FontWeight = 'normal'; % Helvetica字体等转换粗体到Bitmap有问题，直接使用粗体字体名如Helvetica Bold，FontWeight为normal
ax.Title.FontSize = lblfontsize+1;
ax.XLabel.FontName = lblfontname;
ax.YLabel.FontName = lblfontname;
ax.XLabel.FontSize = lblfontsize;
ax.YLabel.FontSize = lblfontsize;
ax.XLabel.Position = [ax.XLabel.Position(1) ax.XLabel.Position(2)+1.08 ax.XLabel.Position(3)];
ax.LineWidth = 0.75;
ax.TickDir = 'none';

% ax2 = axes('Units',get(ax,'Units'),...
%   'Position',get(ax,'Position'),...
%   'Color','none',...
%   'Box','on',...
%   'XColor',"#7E2F8E",...
%   'YColor',"#7E2F8E");
% set(ax2, ...
%   'linewidth',1,...
%   'XTick', [],...
%   'XTickLabel',[],...
%   'YTickLabel',[],...
%   'YTick', []);

svgfilename = ['shapley_' dsname '_ma.svg'];
saveas(fig,svgfilename,'svg');
close(fig);
trimSVGFile(svgfilename);

%% Plot Summary with swarmchart
figure(1); clf; predictorcount = 5;
bars = swarmchart(explainer,NumImportantPredictors=predictorcount,ColorMap='bluered');
fig  = bars(1).Parent.Parent; ax = fig.CurrentAxes;
set(ax,'TickLabelInterpreter','none','FontName','Helvetica Narrow','FontSize',lblfontsize); % global setting for fig.CurrentAxes

yls = yticklabels(ax);
for i=1:length(yls), yls{i} = varnames{str2double(replace(yls{i},'x',''))}; end
yticklabels(ax, yls);

title('');
% ylabel('Predictor (Factor)','FontName',lblfontname,'FontSize',lblfontsize); % ,'Position',[-0.73 predictorcount/2]
ylabel('');
xlabel('Shapely Value','FontName',lblfontname,'FontSize',lblfontsize);
ax.Title.FontName = lblfontname;
ax.Title.FontWeight = 'normal'; % Helvetica字体等转换粗体到Bitmap有问题，直接使用粗体字体名如Helvetica Bold，FontWeight为normal
ax.Title.FontSize = lblfontsize+1;
% ax.YAxis.FontSize = 9; % be setted in global setting for fig.CurrentAxes
% ax.XAxis.FontSize = 9;
ax.XLabel.Position = [ax.XLabel.Position(1) ax.XLabel.Position(2)+1.05 ax.XLabel.Position(3)];
ax.LineWidth = 0.75;
ax.TickDir = 'none';

h = findobj(fig, 'Type', 'colorbar'); 
h.Label.FontName = lblfontname;
h.Label.FontSize = lblfontsize;
h.Label.Position = [h.Label.Position(1)-0.15 h.Label.Position(2) 0];

% svgfilename = ['shapley_' dsname '_sum_swarm.svg'];
% saveas(fig,svgfilename,'svg');
% close(fig);
% trimSVGFile(svgfilename);

%% Plot Summary with boxchart
figure(1); clf;
bars = boxchart(explainer,NumImportantPredictors=5);
fig  = bars(1).Parent.Parent;
ax = fig.CurrentAxes;
set(ax,'TickLabelInterpreter','none','FontName','Helvetica Narrow','FontSize',lblfontsize); % global setting for fig.CurrentAxes

yls = yticklabels(fig.CurrentAxes);
for i=1:length(yls), yls{i} = varnames{str2double(replace(yls{i},'x',''))}; end
yticklabels(fig.CurrentAxes, yls);

ylabel('');
title('');
ax.Title.FontName = lblfontname;
ax.Title.FontWeight = 'normal'; % Helvetica字体等转换粗体到Bitmap有问题，直接使用粗体字体名如Helvetica Bold，FontWeight为normal
ax.Title.FontSize = lblfontsize+1;
ax.XLabel.FontName = lblfontname;
ax.YLabel.FontName = lblfontname;
ax.XLabel.FontSize = lblfontsize;
ax.YLabel.FontSize = lblfontsize;
ax.XLabel.Position = [ax.XLabel.Position(1) ax.XLabel.Position(2)+1 ax.XLabel.Position(3)];
ax.LineWidth = 0.75;
ax.TickDir = 'none';

svgfilename = ['shapley_' dsname '_sum_box.svg'];
saveas(fig,svgfilename,'svg');
close(fig);
trimSVGFile(svgfilename);