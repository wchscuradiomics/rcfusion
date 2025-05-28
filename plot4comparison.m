function fig=plot4comparison(results,legends,titext,precision,style)

if nargin==2,titext=[];precision='%5.4f';style=[]; elseif nargin==3, precision='%5.4f'; style=[]; elseif nargin==4, style=[]; end

if isfield(style,'liwidths')
  liwidths = style.liwidths;
else
  liwidths = 0.5:0.5:length(results)*0.5;
  t = liwidths(fix(length(results)/2)+1:end);
  t = t(length(t):-1:1);
  liwidths(fix(length(results)/2)+1:end) = t;
end

if isfield(style,'colors'), colors=style.colors; else, colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30',...
    '#4DBEEE','#A2142F','r','g','b','c','m','k'}; end
if isfield(style,'specs'), specs=style.specs; else, specs = {'-','-.','--',':','-..'}; end

linestyles = cell(1,length(liwidths));
i = 1; s = 1;
while isempty(linestyles{i})
  linestyles{i} = specs{s};
  if isempty(linestyles{fix(length(results)/2)+i}),  linestyles{fix(length(results)/2)+i} = specs{s}; end
  s = s + 1; i = i + 1;
end
if isempty(linestyles{end}), linestyles{end} = specs{s}; end

aucs = zeros(1,length(results));
for i=1:length(aucs), aucs(i)=results{i}.auc(1); end
[~,aucindices] = sort(aucs);
[liwidths,linewidthindices] = sort(liwidths);
linestyles = linestyles(linewidthindices);

liwidths(aucindices) = liwidths;
linestyles(aucindices) = linestyles;

if length(aucs) > length(colors), error('Please add more colors.'); end
colors = colors(1:length(aucs));
[~,aucdesindices] = sort(aucs,'descend');
colors(aucdesindices) = colors;

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
  plot(results{i}.fprates(:,1),results{i}.tprates(:,1),linestyles{i},'LineWidth', liwidths(i),'Color',colors{i});
  legends{i+1} = [legends{i+1} ' AUC=' num2str(results{i}.auc(1),precision)];
end
pbaspect([0.9999 1 1]);
lg = legend(legends,'Location','southeast');
lg.Position(1) = lg.Position(1) - 0.005;
lg.Position(4) = length(results)*0.067;
xlabel('False Positive Rate (1-Specificity)');
ylabel('True Positive Rate (Sensitivity)');
if ~isempty(titext),title(titext);end
hold off;