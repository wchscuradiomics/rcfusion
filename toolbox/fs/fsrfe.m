function [subset,history] = fsrfe(TRN,trnlabels,cvcvp,nlim,method,seed,parallel)
%[subset,history]=fsrfe(TRN,trnlabels,cvcvp,nlim,method,seed,parallel)

if nargin < 6, seed = 1; end
if nargin < 7, parallel = false; end

rng(seed); % For reproducibility
warning('off');

if ~isequal(unique(trnlabels),[0;1]), error('Only classes [1 0] is supported.'); end

if parallel, maxiter = 100; else, maxiter = 50; end
if strcmpi(method,'corr')
  myfunhandle = @rfecorr;
  opts = statset("Display","off",'UseParallel',parallel);
  rng(1);
  [tf,history] = sequentialfs(myfunhandle,TRN, ...
    "Direction","backward", ...
    "NullModel",true, ...
    "CV","none", ...
    "Options",opts);

  p = size(TRN,2);
  idx = NaN(1,p);
  for i = 1:p, idx(i) = find(history.In(i,:)~=history.In(i+1,:)); end

  threshold = 0.6;
  iter_last_exclude = find(history.Crit(2:end)<threshold,1);
  subset = idx(iter_last_exclude+1:end);
  return;
end

tolfun = 1e-4; tolth = 1e-3;
switch method    
  case 'svm'
    myfunhandle = @rfesvm;
  case 'nb'
    myfunhandle = @rfenb;
  case 'glm'
    myfunhandle = @rfeglm;
    tolfun = 1e-6; tolth = 1e-5; 
  case 'rf'
    myfunhandle = @rferf; 
  case 'gbt'
    myfunhandle = @rfegbt; 
end

opts = statset("Display","off", ...
  "TolFun",tolfun, ...
  "MaxIter",maxiter, ...
  "UseParallel",parallel);
rng(1);
[tf,history] = sequentialfs(myfunhandle,TRN,trnlabels, ...
  "Direction","backward", ...
  "CV",cvcvp, ...
  "Options",opts); 
indices = find(abs(history.Crit-min(history.Crit)) <= tolth);
IN = history.In(indices,:);
crits = history.Crit(indices);
lens = zeros(1,height(IN));
for i=1:length(lens)
  lens(i) = sum(IN(i,:));
end

indices = find(lens == min(lens));
IN = IN(indices,:);
crits = crits(indices);
lens = lens(indices);
[~,index] = min(crits);
subset = find(IN(index(1),:));

if length(subset) > nlim, subset = subset(1:nlim); end
end

function criterion = rfecorr(X)
if size(X,2) < 2
  criterion = 0;
else
  p = size(X,2);
  R = corr(X,"Rows","pairwise");
  R(logical(eye(p))) = NaN;
  criterion = max(abs(R),[],"all");
end
end

function criterion = rfesvm(XTrain,yTrain,XTest,yTest)
  mdl = fitcsvm(XTrain,yTrain,'ClassNames',[1 0]);
  [predictedYTest,score] = predict(mdl,XTest);
  [~,~,~,auc] = perfcurve(yTest,score(:,1),1); criterion = 1 - auc;
  % criterion = sum(~(predictedYTest~=yTest));
end

function criterion = rfenb(XTrain,yTrain,XTest,yTest)
mdl = fitcnb(XTrain,yTrain,'ClassNames',[1 0]);
[predictedYTest,score,~] = predict(mdl,XTest);
[~,~,~,auc] = perfcurve(yTest,score(:,1),1); criterion = 1 - auc;
 % criterion = sum(~(predictedYTest~=yTest));
end

function criterion = rfeglm(XTrain,yTrain,XTest,yTest)
warning('off');
mdl = fitglm(XTrain,yTrain,'Distribution','binomial','Link','logit');
predictedYTest = predict(mdl,XTest);
% [~,~,~,auc] = perfcurve(yTest,predictedYTest(:,1),1);  criterion = 1 - auc;
e = yTest - predictedYTest; criterion = e'*e;
end

function criterion = rferf(XTrain,yTrain,XTest,yTest)
mdl = TreeBagger(100, XTrain, yTrain, ...
  'Method', 'classification', ...
  'OOBPrediction', 'On');
[predictedYTest,score] = predict(mdl,XTest);
[~,~,~,auc] = perfcurve(yTest,score(:,1),1); criterion = 1 - auc;
 % criterion = sum(~(predictedYTest~=yTest));
end

function criterion = rfegbt(XTrain,yTrain,XTest,yTest)
t = templateTree('MaxNumSplits',5);
mdl = fitcensemble(XTrain,yTrain, ...
  'Method','AdaBoostM1', ...
  'Learners',t, ...
  'ClassNames',[1 0]);
[predictedYTest,score] = predict(mdl,XTest);
[~,~,~,auc] = perfcurve(yTest,score(:,1),1); criterion = 1 - auc;
 % criterion = sum(~(predictedYTest~=yTest));
end