clear; clc;
fsnames = {'swr.mat','stdga.mat','relieff.mat','mrmr.mat','lasso.mat','inf.mat','il.mat','rfe.mat','our.mat'}; 
format = '%4.3f';

for i=1:length(fsnames)
  fsname = fsnames{i};

  param = getfield(load(fsname,'param'),'param');
  li = getfield(load(fsname,'li'),'li');
  LCL = getfield(load(fsname,'LCL'),'LCL');
  LCLSCNO4H = getfield(load(fsname,'LCLSCNO4H'),'LCLSCNO4H');
  tstresult = getfield(load(fsname,'tstresult'),'tstresult');

  radscores = LCL(li.tstindices,:)*param.fs.lcl.coefs + param.fs.lcl.intercept;
  scno4hradscores = LCLSCNO4H*param.fs.lcl.coefs + param.fs.lcl.intercept;
  radresult = perfresult(li.tstlabels,radscores,param.performance.version);
  scno4hradresult = perfresult(li.scno4hlabels,scno4hradscores,param.performance.version);

  disp([ ...
    num2str(radresult.auc(1),format) ' (' num2str(radresult.auc(2),format) '-' num2str(radresult.auc(3),format) ')' 9 ...
    num2str(radresult.accuracy(1),format) ' (' num2str(radresult.accuracy(2),format) '-' num2str(radresult.accuracy(3),format) ')' 9 ...
    num2str(scno4hradresult.auc(1),format) ' (' num2str(scno4hradresult.auc(2),format) '-' num2str(scno4hradresult.auc(3),format) ')' 9 ...
    num2str(scno4hradresult.accuracy(1),format) ' (' num2str(scno4hradresult.accuracy(2),format) '-' num2str(scno4hradresult.accuracy(3),format) ')' ...
  ]);

  % disp(num2str((radresult.auc(1)*(128/188)+scno4hradresult.auc(1)*(60/188))));
end