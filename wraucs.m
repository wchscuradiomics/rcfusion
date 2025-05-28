function wraucs(methods,names,li)

form = '%4.3f';

for i=1:length(methods)
  s = names{i};
  valresult = getfield(load(methods{i},'valresult'),'valresult','nom');
  valresult.labels = li.trnlabels;
  tstresult = getfield(load(methods{i},'tstresult'),'tstresult','nom');
  tstresult.labels = li.tstlabels;
  scno4hresult = getfield(load(methods{i},'scno4hresult'),'scno4hresult','nom');
  scno4hresult.labels = li.scno4hlabels;
  results = {valresult,tstresult,scno4hresult};
  for j=1:length(results)
    r = results{j};
    r = perfresult(r.labels,r.SCORE(:,1),'optimal:OPTROCPT|nboot:100');
    s = [s 9 num2str(r.auc(1),form) ' (' num2str(r.auc(2),form) ' - ' ...
      num2str(r.auc(3), form) ')' 9 ... 
    num2str(r.accuracy(1),form) ' (' num2str(r.accuracy(2),form) ' - ' ...
    num2str(r.accuracy(3), form) ')' 9]; %#ok<AGROW> 
  end  
  disp(s);
end