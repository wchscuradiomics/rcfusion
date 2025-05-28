function subset = mdl2subset(mdl)
subset = zeros(1,length(mdl.PredictorNames));
for i=1:length(subset)
  subset(i) = str2double(replace(mdl.PredictorNames{i},'x',''));
end
end