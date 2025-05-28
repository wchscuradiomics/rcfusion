function subset = rmvrdnfeats3pearsoncc(TRN,th)
%subset=rmvrdnfeatures(TRN) removes redundant features using Pearson's correlation coefficients. It is applicable to situations
%where there is a linear relationship between two continuous variables and the data are required to be approximately normally
%distributed (when the sample size is large, it can be considered to be approximately normally distributed).
%
% TRN: a n*m matrix specifying n samples with m features (continuous variables are preferred).
%
% th: a double representing the threshold. One of the two features will be removed if the correlation coefficient between the two
% features is > th.

if nargin == 1, th = 0.8; end % the default threshold

R = abs(corrcoef(TRN));
avgrs = mean(R);
[a,b] = meshgrid(1:size(R,1),1:size(R,1));
R(a<=b) = nan;

rmvindices = zeros(length(avgrs)^2,1);
k = 0;
for i=1:size(R,1)
  for j=1:size(R,2)
    if isnan(R(i,j)), continue; end
    if R(i,j) > th
      k = k + 1;
      if avgrs(i) > avgrs(j),rmvindices(k) = i; else, rmvindices(k) = j; end
    end
  end
end
rmvindices = rmvindices(1:k);
rmvindices = unique(rmvindices);

subset = 1:size(TRN,2);
subset(rmvindices) = [];