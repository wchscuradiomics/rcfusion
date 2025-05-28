function n = auton1scores(sorscores,method,nlim)
%n=auton1scores(sorscores,method,nlim) automatically select n based on ranked scores.

if ~issorted(sorscores,'descend'), error('sorscores must be a descend sorted vector.'); end
if nlim < 3 || length(sorscores) < 3, error('nlim or length(sorscores) must be >= 5.'); end

if contains(method,'drop','IgnoreCase',true) % recommended
  if strcmpi(method,'drop'), s = 3; else, s = str2double(replace(method,'drop','')); end
  if s > nlim, s = nlim; end

  dropheights = abs(diff([sorscores sorscores(end)])); 
  smooths = abs(diff([sorscores sorscores(end) sorscores(end)],2)); 
  dropprecision = min([10^( floor(log10(abs(sorscores(1)))) )/10 0.001]);
  if abs(sorscores(1)) < 10
    smthprecision = 10^( floor(log10(abs(sorscores(1)))) )/10;
  else
    smthprecision = 10^( floor(log10(abs(sorscores(1))))-1 )/10;
  end
  
  cumscores = cumsum(sorscores/sum(sorscores));
  st = max([find(cumscores >= 0.5,1,'first') 3 s]);

  [~, sordropindices] = sort(dropheights,'descend');
  maxdropindex = sordropindices(find(sordropindices >= st,1,"first"));

  if maxdropindex >= nlim,n = maxdropindex; return; end

  for i=maxdropindex:length(sorscores)-1
    if dropheights(i+1) <= dropprecision && smooths(i+1) <= smthprecision, n = i; break; end
  end
  if n > nlim, n = nlim; end

elseif contains(method,'top','IgnoreCase',true)
  if strcmpi(method,'top'), n = nlim; else, n = str2double(replace(method,'top','')); end
  if n > nlim, n = nlim; end

elseif contains(method,'cumscore','IgnoreCase',true) % be used with caution
  th = str2double(replace(method,'cumscore',''));
  cumscores = cumsum(sorscores/sum(sorscores));
  n = find(cumscores >= th,1,'first');

  if n > nlim, n = nlim; end
else
  error(['Not supported method: ' method]);
end

