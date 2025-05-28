function [predicted,observed,pvalue,chi2] = hl(probabilities,labels,gpcount) % Hosmer–Lemeshow test

[probabilities,sortindices] = sort(probabilities);
labels = labels(sortindices);

observed = zeros(1,gpcount);
predicted = zeros(1,gpcount);
chi2 = 0;
percentiles = [-1 prctile(probabilities, (1/gpcount:1/gpcount:1)*100)]; 
indicess = cell(gpcount,1);
for i=1:gpcount
  indicess{i} = find(probabilities > percentiles(i) & probabilities <= percentiles(i+1));
end

for i=1:gpcount
  indices = indicess{i};

  smpcount = length(indices); % the number of elements in the i-th group: count of
  obspossmpcount = sum(labels(indices)==1);

  observed(i) = obspossmpcount/smpcount;
  predicted(i) = mean(probabilities(indices));

  predpossampcount = sum(probabilities(indices)); % predicted(i) * smpcount; sum(probabilities(indices))
  chi2 = chi2 + (obspossmpcount - predpossampcount)^2/predpossampcount+ ...
     ((smpcount-obspossmpcount) - (smpcount-predpossampcount))^2/(smpcount-predpossampcount);
end

df = gpcount - 2;
pvalue = 1 - chi2cdf(chi2, df);