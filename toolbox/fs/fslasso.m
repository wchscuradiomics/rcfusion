function [subset,mse,B,fitinfo] = fslasso(TRN,trnlabels,cvcvp,preoption,standardize,criteria,alpha, paroption,dfmax)
%[subset,mse,B,fitinfo] = fslasso(TRN,trnlabels,cvcvp,preoption,standardize,criteria,alpha, paroption,dfmax) perform feature
%selection using LASSO.
%
% A row vector of TRN is a sample.
%
% preoption: refer to fs3multistages. For example, 'utest|corrcoef-remove|0.05|0.8', 'none|nan', ....
%
% criteria: 'MSE' | '1SE'.

%% check and set parameters
if nargin < 8, preoption = false; end
if nargin < 9, dfmax = Inf;  end
if isa(cvcvp,'double'), cvcvp = cvpartition(trnlabels,'KFold',cvcvp); end

strs = strsplit(preoption,'|');
if contains(preoption,'two-stages') % old version
  premethods = cell(1,2);
  if contains(strs{1},'utest')
    premethods{1} = 'utest';
  else
    premethods{1} = 'ttest';
  end
  premethods{2} = 'corrcoef-remove';
  ths = str2double(strs(2:3));  
else
  premethods = strs(1:length(strs)/2);
  ths = str2double(strs(length(strs)/2 + 1:end));
end
%% preselection
presubset = fs3multistages(TRN,trnlabels,premethods,ths);
TRNSUBSET = TRN(:,presubset);

%% LASSO
[B,fitinfo] = classo(TRNSUBSET,trnlabels, ...
  'CV',cvcvp, ...
  'Alpha' ,alpha, ...
  'Options',statset('UseParallel',paroption),...
  'Standardize',standardize,'DFmax',dfmax);
if strcmpi(criteria,'none')
  subset = nan; mse = nan;
else
  fitinfo.bstlindex = selbstlambda(fitinfo, criteria);
  subset = presubset(abs(B(:,fitinfo.bstlindex)) > 1e-10);
  mse = fitinfo.MSE(fitinfo.bstlindex);
end
