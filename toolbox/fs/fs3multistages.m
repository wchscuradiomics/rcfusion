function subset = fs3multistages(TRN,tralabels,methods,ths) 
%subset=fs3multistages(TRN,tralabels,methods,ths) performs feature selection using multiple stages.
%
% TRN: a n-by-m matrix representing training samples, where n is the number of training samples and m is the number of features.
%
% tralabels: a n-by-1 vector representing the lables of TRN.
%
% methods: a 1-by-c cell representing feature selection methods, where c is the number of methods.
%
% ths: a 1-by-c vector representing the threshold values of methods.
%
% For example:
% presubset = fs3multistages(TRN,tralabels,{'ttest','corrcoef-remove'},[0.05 0.8]). Feature selection methods refer to function
% fs31stage (type help fs31stage in a MATLAB command window).

subsets = cell(1,length(ths));
for c = 1:length(ths)  
  subsets{c} = fs31stage(TRN,tralabels,methods{c},ths(c));  
  TRN = TRN(:,subsets{c});
end

if length(subsets) == 1
  subset = subsets{1};
else
  subset = subsets{end};
  for c = length(subsets):-1:2    
    cnsub1 = subsets{c-1};
    subset = cnsub1(subset);
  end
end

end % end for function fs3multistages