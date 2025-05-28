function [subset,idx,scores] = fsrelieff(varargin)
%[subset,idx,scores]=fsrelieff(TRN,trnlabels,k,method,way,nlim,verbose)
%
% k: the number of neighbors.
%
% method: 'topx'|'drop'.

[TRN,trnlabels,k,method,way] = varargin{1:5};
if nargin >= 6, nlim = varargin{6}; else, nlim = fix(height(TRN)/10); end
if nargin >= 7, verbose = varargin{7}; else, verbose = false; end

if nlim>width(TRN), nlim=width(TRN); end

[idx,scores] = relieff(TRN, trnlabels, k, 'method',method);
sorscores = scores(idx);
if verbose, bar(sorscores); end

% The drop in score between the i-th and (i+1)-th most important predictors is large, while the drops after the j (j > i+1)
% predictor are relatively small. A drop in the importance score represents the confidence of feature selection. Therefore, the
% large drop implies that the software is confident of selecting the most important predictor. The small drops indicate that the
% difference in predictor importance are not significant.
dropindex = auton1scores(sorscores,way,nlim);
if isempty(dropindex)
  subset = idx(1:nlim);
elseif dropindex > nlim
	if dropindex <= fix(height(TRN)/5), subset = idx(1:dropindex); else, subset = []; end
else
  subset = idx(1:dropindex);
end % end for function fsrelieff