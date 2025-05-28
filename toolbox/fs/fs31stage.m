function subset = fs31stage(D,labels,method,th) 
%subset=fs31stage(D,labels,method,th) selects features with one stage (ttest, utest, Pearson's correlation coefficients for
%removing redundant features, and Spearman's correlation coefficient for outcome/label, univariate regression of p-value,
%univariate regression of R-squared), or variable coefficient (cv).
%
% method: 'ttest' | 'utest' | 'corrcoef-remove' | 'corrcoef-label' | 'ur-pvalue' | 'ur-rsquared' | 'cv' | 'none', a string
% specifying an abbreviative method name
%
% th: a double number representing the threshold

names = unique(labels);
if length(names) ~=2, error('There must be two types in labels.'); end

switch method
  case 'ttest'
    ps = nan(1,size(D,2));
    if isa(D,'table')
      for j=1:size(D,2)
         v = D.(j);
         if isa(v,'categorical')
           tt = zeros(size(labels));
           tt(labels==names(1)) = 1;
           tt(labels==names(2)) = 2;
           tt = categorical(tt);           
           [ps(j),~] = fisherorchi(tt,v);
         end 
      end % end for TRN is a table
    else % TRN is a double matrix
      for j=1:size(D,2)
        v = D(:,j);
        [~,ps(j)] = ttest2(v(labels==names(1)),v(labels==names(2)));        
      end
    end
    subset = find(ps<=th);
  case 'utest'
    ps = nan(1,size(D,2));
    if isa(D,'table')
      for j=1:size(D,2)
         v = D.(j);
         if isa(v,'categorical')
           tt = zeros(size(labels));
           tt(labels==names(1)) = 1;
           tt(labels==names(2)) = 2;
           tt = categorical(tt);           
           ps(j) = fisherorchi(tt,v);
         end 
      end % end for TRN is a table
    else % TRN is a matrix
      for j=1:size(D,2)
        v = D(:,j);
        [ps(j),~,~] = ranksum(v(labels==names(1)),v(labels==names(2)));
      end      
    end
    subset = find(ps<=th);
  case 'corrcoef-remove'
    subset = rmvrdnfeats3pearsoncc(D,th);
  case 'corrcoef-label' % Spearman's correlation coefficient
    pvals = zeros(1,size(D,2));
    for j=1:size(D,2)
      [~,pvals(j)] = corr(D(:,j),labels,'Type','Spearman');
    end
    subset = find(pvals <= th);
  case 'ur-pvalue'
    ps = nan(1,size(D,2));
    for j=1:size(D,2)
      % mdl =  fitglm(TRN(:,j),trnlabels,'Distribution','binomial'); % fitlm(TRN(:,j),trnlabels);
      v = D(:,j);
      v = addvars(v,labels);
      mdl = fitglm(v,'linear','Distribution','binomial','Link','logit');
      [ps(j),~] = coefTest(mdl);
    end
    subset = find(ps <= th);
  case 'ur-rsquared'
    rs = nan(1,size(D,2));
    for j=1:size(D,2)
      mdl =  fitglm(D(:,j),labels,'Distribution','binomial');
      [rs(j),~] = mdl.Rsquared.Ordinary;
    end
    subset = find(rs >= th);
  case 'cv'
    D = D - min(D);
    v = std(D)./ mean(D);
    subset = find(v < th);
  case 'none'
    subset = 1:size(D,2);
  otherwise
    error([method ' is not implemented now.']);
end

end % end for function fs31stage