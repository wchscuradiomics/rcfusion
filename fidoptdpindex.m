function dpindex = fidoptdpindex(sorscores,smooth,nlim)
%dpindex=fidoptdpindex(sorscores,nlim,smooth) finds optimized drop index.

if nlim < length(sorscores), sorscores = sorscores(1:nlim); end

dropheights = abs(diff(sorscores)); % length(dropheights) equals to length(sorscores)-1
dropheights = dropheights(1:end-2);

copyscores = [sorscores sorscores(end) sorscores(end)];
d2s = zeros(1,length(dropheights)+1);
for j=1:length(dropheights)+1
  d2s(j) = copyscores(j+2) - 2*copyscores(j+1) + copyscores(j);
end
d2s = abs(d2s(2:end));

if length(d2s) ~= length(dropheights), error('No consistent between d2s and dropheights.'); end

candindices = find(d2s <= smooth);

[~, candindex] = max(dropheights(candindices) - d2s(candindices));
dpindex = candindices(candindex);

end