%% Train
clear; clc;
mlnames = {'svm','nb','rf','lgb','xgb'}; % 
fsnames = {'swr','stdga','relieff','mrmr','lasso','inf','il','rfe','our'}; % 

for i = 1:length(fsnames)
  % if exist([fsnames{i} '_ml.mat'],'file'), continue; end
  li = getfield(load([fsnames{i} '.mat'],'li'),'li');
  TRN = getfield(load([fsnames{i} '.mat'],'TRN'),'TRN');
  TST = getfield(load([fsnames{i} '.mat'],'TST'),'TST');
  TSTSCNO4H = getfield(load([fsnames{i} '.mat'],'TSTSCNO4H'),'TSTSCNO4H');
  param = setclparam(mlnames,struct('suffix',[fsnames{i} '_ml'],'varcount',width(TRN)));
  [classifier,trnresult,valresult,~,~] = convclassify(TRN,li,param,[fsnames{i} '_ml'], mlnames);
  for j = 1:length(mlnames)
    tstscores = cpredict(classifier.(mlnames{j}),TST);
    tstscoresscno4h = cpredict(classifier.(mlnames{j}),TSTSCNO4H);
    tstresult.(mlnames{j}) = perfresult(li.tstlabels,tstscores,'optimal:youden|nboot:100');
    scno4hresult.(mlnames{j}) = perfresult(li.scno4hlabels,tstscoresscno4h,'optimal:youden|nboot:100');
  end

  save([fsnames{i} '_ml.mat'], ...
    'li','TRN','TST','TSTSCNO4H','classifier',...
    'trnresult','valresult','tstresult','scno4hresult','param' ...
  );
end

%% Print
clear; clc;
fsnames = {'swr','stdga','relieff','mrmr','lasso','inf','il','rfe','our'};
mlnames = {'svm','nb','rf','lgb','xgb'};
format = '%4.3f';
for i = 1:length(fsnames)
  trnresult = getfield(load([fsnames{i} '_ml.mat'],'trnresult'),'trnresult');
  valresult = getfield(load([fsnames{i} '_ml.mat'],'valresult'),'valresult');
  tstresult = getfield(load([fsnames{i} '_ml.mat'],'tstresult'),'tstresult');
  scno4hresult = getfield(load([fsnames{i} '_ml.mat'],'scno4hresult'),'scno4hresult');
  disp(upper(fsnames{i}));
  for j=1:length(mlnames)
    s = [upper(mlnames{j}) 9 ...
      num2str(trnresult.(mlnames{j}).auc(1),format) ' (' ...
      num2str(trnresult.(mlnames{j}).auc(2),format) '-' num2str(trnresult.(mlnames{j}).auc(3),format) ')' 9 ...
      num2str(trnresult.(mlnames{j}).accuracy(1),format) ' (' ...
      num2str(trnresult.(mlnames{j}).accuracy(2),format) '-' num2str(trnresult.(mlnames{j}).accuracy(3),format) ')' 9 ...
      num2str(valresult.(mlnames{j}).auc(1),format) ' (' ...
      num2str(valresult.(mlnames{j}).auc(2),format) '-' num2str(valresult.(mlnames{j}).auc(3),format) ')' 9 ...
      num2str(valresult.(mlnames{j}).accuracy(1),format) ' (' ...
      num2str(valresult.(mlnames{j}).accuracy(2),format) '-' num2str(valresult.(mlnames{j}).accuracy(3),format) ')' 9 ...
      num2str(tstresult.(mlnames{j}).auc(1),format) ' (' ...
      num2str(tstresult.(mlnames{j}).auc(2),format) '-' num2str(tstresult.(mlnames{j}).auc(3),format) ')' 9 ...
      num2str(tstresult.(mlnames{j}).accuracy(1),format) ' (' ...
      num2str(tstresult.(mlnames{j}).accuracy(2),format) '-' num2str(tstresult.(mlnames{j}).accuracy(3),format) ')' 9 ...
      num2str(scno4hresult.(mlnames{j}).auc(1),format) ' (' ...
      num2str(scno4hresult.(mlnames{j}).auc(2),format) '-' num2str(scno4hresult.(mlnames{j}).auc(3),format) ')' 9 ...
      num2str(scno4hresult.(mlnames{j}).accuracy(1),format) ' (' ...
      num2str(scno4hresult.(mlnames{j}).accuracy(2),format) '-' num2str(scno4hresult.(mlnames{j}).accuracy(3),format) ')'];
    disp(s);
  end
  disp('');
end

%%
clear;clc;
mlnames = {'svm','nb','rf','lgb','xgb'};
fsnames = {'swr','stdga','relieff','mrmr','lasso','inf','il','rfe','our'};
for i=1:length(fsnames)
  clearvars -except mlnames fsnames i j;
  load([fsnames{i} '_ml.mat']);
  for j = 1:length(mlnames)
    if strcmpi(mlnames{j},'lgb') || strcmpi(mlnames{j},'xgb')
      classifier.(mlnames{j}) = py.rwdata.pickleread(param.(mlnames{j}).varpath);
    end
    tstscores = cpredict(classifier.(mlnames{j}),TST);
    tstresult.(mlnames{j}) = perfresult(li.tstlabels,tstscores,'optimal:youden|nboot:100');
  end

  save([fsnames{i} '_ml.mat'], ...
    'li','TRN','TST','TSTSCNO4H','classifier',...
    'trnresult','valresult','tstresult','scno4hresult','param' ...
  );
end