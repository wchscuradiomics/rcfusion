function [pvalues,rsvalues,subset] = univaranalysis(T,labels,method,threshold)
%[pvalues,rsvalues,subset]=univaranalysis(T,labels,method,threshold) select features using univariate analysis.
%
% T: a n-by-m table or matrix with n samples and m variables. If categorical variables exist, please be sure correct types are
% already in T for these categorical variables.
%
% labels: a n-by-1 representing 0 or 1 for n samples.
%
% method: 'regression' | 'diff' (significant difference test)
%
% threshold: a double scalar in (0, 1) to select a subset (features).
%
% pvalues: a 1-by-m double vector representing p-values of m variables.
%
% rsvalues: a 1-by-m double vector representing R-square of m variables in regression. If 'diff' is used, rsvalues is nan.
%
% subset: a 1-by-x double vector repsenting x selected variables based on threshold and method arguments.

pvalues = nan(1,width(T));
rsvalues = nan(1,width(T));
if strcmpi(method,'regression')
  for j=1:size(T,2)
    if isa(T,'table')
      mdl = fitglm(addvars(T(:,j),labels),'linear', ...
        'Distribution','binomial', ...
        'Link','logit');
    else
      mdl = fitglm(T(:,j),labels,'linear', ...
        'Distribution','binomial', ...
        'Link','logit');
    end    
    [pvalues(j),~] = coefTest(mdl);
    rsvalues(j) = sqrt(mdl.Rsquared.Adjusted);
  end
  subset = find(pvalues <= threshold);
else % significant difference test  
  type = 'chi'; % type = fisherorchi(T, 2);
  clnames = unique(labels);
  for j=1:size(T,2)
    if isa(T,'table'), v = T.(j); else, v=T(:,j); end
    if isa(v,'categorical')
      tt = NaN(length(labels),1);
      tt(labels==clnames(1)) = 1;
      tt(labels==clnames(2)) = 2;
      tt = categorical(tt);      
      if strcmpi(type,'fisher')
        [~,pvalues(j),~] =  fishertest(crosstab(tt,v));        
      else % chi
        [~,~,pvalues(j)] = crosstab(tt,v);
      end
    else % v is a double vector
      pvalues(j) = ranksum(v(labels==clnames(1)),v(labels==clnames(2)));
      % [~,pvalues(j)]=ttest2(v(labels==clnames(1)),v(labels==clnames(2)));
    end    
  end
  subset = find(pvalues <= threshold);
end