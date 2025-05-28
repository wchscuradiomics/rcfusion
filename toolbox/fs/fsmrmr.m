function [subset,idx,scores] = fsmrmr(TRN,trnlabels,preoption,way,nlim,verbose)
%[subset,idx,scores]=fsmrmr(TRN,trnlabels,preoption,way,nlim,verbose) perform feature selection using MRMR.
%
% A row vector of TRN is a sample.
%
% preoption: refer to fs3multistages. For example, 'utest|corrcoef-remove|0.05|0.8', 'none|nan', ....
%
% fscmrmr: a matlab function that ranks features for classification using the minimum redundancy maximum relevance (MRMR)
% algorithm.
%
% Use the drop in score between the i-th and (i+1)-th most important predictors --> select n features, where n <= nlim.
%
% Note: if drop or topn ways are not effective, try a topi~topn way.

%% check and set parameters
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

if nargin < 5, nlim = fix(size(TRN,1)/10); end
if nargin < 6, verbose = true; end
if nlim > width(TRN), nlim = width(TRN); end

%% preselection & set parameters of nlim and difference
presubset = fs3multistages(TRN,trnlabels,premethods,ths);
TRNSUBSET = TRN(:,presubset);
if nlim > size(TRNSUBSET,2), nlim = size(TRNSUBSET,2); end

%% MRMR
[idx,scores] = fscmrmr(TRNSUBSET, trnlabels);
sorscores = scores(idx);
if verbose, bar(sorscores); end

% The drop in score between the i-th and (i+1)-th most important predictors is large, while the drops after the j (j > i+1)
% predictor are relatively small. A drop in the importance score represents the confidence of feature selection. Therefore, the
% large drop implies that the software is confident of selecting the most important predictor. The small drops indicate that the
% difference in predictor importance are not significant.
dropindex = auton1scores(sorscores,way,nlim);
if isempty(dropindex)
  subset = presubset(idx(1:nlim));
elseif dropindex > nlim
  if dropindex <= fix(height(TRN)/5), subset = presubset(idx(1:dropindex)); else, subset = []; end
else
  subset = presubset(idx(1:dropindex));
end % end for function fsmrmr