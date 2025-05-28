function [subset,idx,scores] = fsrelieff(TRN,trnlabels,k,method,way,nlim,verbose)
%[subset,idx,scores]=fsrelieff(TRN,trnlabels,k,method,way,nlim,verbose)
%
% k: the number of neighbors.
%
% method: 'topx'|'drop'.

[idx,scores] = relieff(TRN, trnlabels, k, 'method',method);
sorscores = scores(idx);
if verbose, bar(sorscores); end

if contains(way,'top') % get the top n features
  topn = str2double(replace(way,'top',''));
  if length(idx)>=topn
    subset = idx(1:topn);
  else
    subset = idx(1:end);
  end
  return;
end

% The drop in score between the i-th and (i+1)-th most important predictors is large, while the drops after the j (j > i+1)
% predictor are relatively small. A drop in the importance score represents the confidence of feature selection. Therefore, the
% large drop implies that the software is confident of selecting the most important predictor. The small drops indicate that the
% difference in predictor importance are not significant.
dropindex=fidrop(sorscores,nlim);

if isempty(dropindex)
  subset = idx(1:nlim);
else
  subset = idx(1:dropindex);
end % end for function fsrelieff