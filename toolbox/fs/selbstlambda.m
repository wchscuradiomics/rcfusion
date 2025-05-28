function bstlindex = selbstlambda(fitinfo,criteria)
%bstlindex=selbstlambda(fitinfo,criteria) selects the best lambda index.

criteria = char(criteria);

if strcmpi(criteria,'MSE')
  bstlindex = fitinfo.IndexMinMSE;
elseif strcmpi(criteria,'CVAUC')
  [~, bstlindex] = max(fitinfo.cvaucs);
elseif strcmpi(criteria,'MAUC')
  bstlindex = fitinfo.IndexMaxMAUC;
elseif contains(criteria,'MAUC') % IndexSEMAUCs
  i = str2double(criteria(5));
  if i > 3
    error('<=3 SE MAUCs are allowed.');
  elseif i == 1
    bstlindex = fitinfo.Index1SEMAUC;
  else
    bstlindex = fitinfo.IndexSEMAUCs(i);
  end
elseif contains(criteria,'MSE') % IndexSEs
  i = str2double(criteria(4));
  if i > 3
    error('<=3 SEs are allowed.');
  elseif i == 1
    bstlindex = fitinfo.Index1SE;
  else
    bstlindex = fitinfo.IndexSEs(i);
  end
else
  error('Invalid criteria.');
end