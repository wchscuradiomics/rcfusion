function fig = plotcalcurves(varargin)
%plotcalcurves(predicted,observed,legnames, ... xlabstring,ylabstring,titstring) plots a calibration curve.
%
% For example, 
% titstring = 'Calibration Curve for Model(s) XXX'; 
% xlabstring = 'Estimated Probability';
% ylabstring = 'Observed Probability';
% legnames = {'Idea','SVM','RF'};

[predicted,observed,legnames,xlabstring,ylabstring,titstring] = varargin{:};

if ~isa(predicted,'cell') || ~isa(observed,'cell') || ~isa(legnames,'cell')
  error('Params predicted, observed, and legend must be cell.');
end

markers = {'*','+','o','s','d','^','v','>','<','p','h','_','|',};
colors = {'#0072BD','#D95319','#EDB120','#7E2F8E','#77AC30','#4DBEEE', '#A2142F','r','g','b','c','m','y'};
plot([0 1], [0 1], 'LineStyle','--','Color','#999999','LineWidth',0.75);
for i=1:length(predicted)
  hold on;
  if i~=length(predicted)
    plot(predicted{i},observed{i},'Color',colors{i},'Marker',markers{i}, ...
      'LineStyle','-','LineWidth',0.75); 
  else % '*',
    plot(predicted{i},observed{i},'Color','#800080','Marker','x', ...
      'LineStyle','-','LineWidth',1.5);
  end
end
legend(legnames,'Location',[0.676 0.176 0.2 0.046*(length(legnames)+1)]);
title(titstring);
xlabel(xlabstring);
ylabel(ylabstring);
hold off;

legend boxoff;