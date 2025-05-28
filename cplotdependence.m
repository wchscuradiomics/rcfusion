function fig = cplotdependence(explainer,varindex,clrvarindex,varnames,lblfontsize,lblfontname)

varname = varnames{varindex};

% lblfontsize = 9; lblfontname = 'Helvetica';

sca = plotDependence(explainer,varindex,ColorPredictor=clrvarindex);
fig = sca.Parent.Parent;
ax = fig.CurrentAxes;
set(ax,'TickLabelInterpreter','none','FontName','Helvetica Narrow','FontSize',lblfontsize); % global setting for fig.CurrentAxes

ylabel(ax,['Shapley Values for ' varname],'FontName',lblfontname,'FontSize',lblfontsize);
xlabel(ax,varname,'FontName',lblfontname,'FontSize',lblfontsize);
title(ax,['Dependence Plot for ' varname ' and ' varnames{clrvarindex}], ...
  'FontName',lblfontname,'FontWeight','normal','FontSize',lblfontsize+1);
title(ax,'');

h = findobj(fig, 'Type', 'colorbar'); 
h.Label.String = varnames{clrvarindex};
h.Label.FontName = lblfontname;
h.Label.FontSize = lblfontsize;
h.Label.Position = [0 h.Label.Position(2) h.Label.Position(3)];

ax.YLabel.Position = [ax.YLabel.Position(1)+1.3 ax.YLabel.Position(2) ax.YLabel.Position(3)];
ax.XLabel.Position = [ax.XLabel.Position(1) ax.XLabel.Position(2)+0.25 ax.YLabel.Position(3)];
ax.LineWidth = 0.75;
ax.TickDir = 'none';