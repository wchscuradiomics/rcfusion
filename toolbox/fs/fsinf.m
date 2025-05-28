function [subset,idx,scores] = fsinf(varargin)
%[subset,idx,scores]=fsinf(TRN,trnlabels,alphas,us,preoption,nlim,verbose) performs feature selection using InfFS_U or InfFS_S.
%
% preoption: 'two-stages' | 'correlation' | 'ttest' or its extension that includes threshold parameters of p-value and correlation
% value, i.e., 'xxx | pthreshold | cthreshold'. For example, 'two-stages' or 'two-stages|0.05|0.8'.

% check and set parameters
[TRN,trnlabels,alphas,us,preoption] = varargin{1:5};
strs = strsplit(preoption,'|');
premethods = strs(1:length(strs)/2);
ths = str2double(strs(length(strs)/2 + 1:end));
if nargin >= 6, nlim = varargin{6}; else, nlim = fix(size(TRN,1)/10); end
if nargin >= 7, verbose = varargin{7}; else, verbose = true; end

% preselection & set parameters of nlim and difference
presubset = fs3multistages(TRN,trnlabels,premethods,ths);
TRNSUBSET = TRN(:,presubset);
if nlim > size(TRNSUBSET,2), nlim = size(TRNSUBSET,2); end

% InfFS
if strcmpi(us, 'U')
  [idx, scores, subindices] = InfFS_U(TRNSUBSET, trnlabels, alphas(1), verbose);
elseif strcmpi(us, 'S')
  [idx, scores, subindices] = InfFS_S(TRNSUBSET, trnlabels, alphas);
else
  error(['Not supported way: ' prelexptparam.way]);
end

subset = presubset(subindices);
if nlim < length(subset), subset=subset(1:nlim); end

if verbose, bar(sorscores); end