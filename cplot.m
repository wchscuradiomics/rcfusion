function fig = cplot(results,legends,titext,precision)

if nargin == 2
  titext = []; precision='%5.4f';
elseif nargin==3
  precision = '%5.4f';
end

if length(results) > 5
  liwidths = 0.25:0.25:0.25*length(results);
else
  liwidths = 0.25:0.5:(0.5*length(results)+0.25);
end

colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE', ...
  '#A2142F'};

fig = figure;
hold on;
plot(0:0.01:1,0:0.01:1,'--','Color','#dddddd');
box on;
legends = [{''} legends];
for i=1:length(results)
  [m, n] = size(results{i}.fprates);
  if m < n
   results{i}.fprates = results{i}.fprates';
   results{i}.tprates = results{i}.tprates';
  end
  hold on;
  plot(results{i}.fprates(:,1),results{i}.tprates(:,1), ...
    'LineWidth', liwidths(i), 'Color',colors{i});
  legends{i+1} = [legends{i+1} ' AUC=' ...
    num2str(results{i}.auc(1),precision)];
end
pbaspect([0.9999 1 1]);
lg = legend(legends,'Location','southeast');
lg.Position(4) = length(results)*0.1;
xlh = xlabel('False Positive Rate (1 - Specificity)');
ylh = ylabel('True Positive Rate (Sensitivity)');
if ~isempty(titext),tih=title(titext); else, tih=[]; end
hold off;

fontname = 'Helvetica';
ax=fig.CurrentAxes;
set(ax,'TickLabelInterpreter','none','FontName',fontname,'FontSize',10);
ax.YAxis.FontSize = 9;
ax.XAxis.FontSize = 9;
set(gca,'xTick',0:0.1:1);
set(gca,'xTickLabel',num2str(get(gca,'xTick')','%.1f'));
set(gca,'yTick',0.1:0.1:1);
set(gca,'yTickLabel',num2str(get(gca,'yTick')','%.1f'));
lg.FontSize = 10;
lg.FontName = fontname;
xlh.FontSize = 10;
xlh.FontName = fontname;
ylh.FontSize = 10;
ylh.FontName = fontname;
xlh.Position = [0.5 -0.03];
ylh.Position = [-0.06 0.5];
if ~isempty(tih)
  tih.FontSize = 12;
  tih.FontName = fontname;
  % tih.FontWeight = 'normal';
end
a = gca;
 % negative numbers move the ticklabels down (positive -> up)
a.XRuler.TickLabelGapOffset = -3;   
% negative numbers move the ticklabels right (negative -> left)
a.YRuler.TickLabelGapOffset = -1;    