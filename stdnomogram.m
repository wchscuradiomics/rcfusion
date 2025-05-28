function [x,util] = stdnomogram(CAND,CLN,li,param)

% subset based on ga
rng(param.seed);
[x,util] = stdga(CAND(li.trnindices,:),CLN(li.trnindices,:),li.trnlabels,li.cvcvp,param); 
isxdisplayed = false;
for iindividual = 1:height(util.X)
  if (isequal(util.X(iindividual,:),x) && isxdisplayed) || ~isequal(util.X(iindividual,:),x),continue; end
  % if (isequal(util.X(iindividual,:),x) && isxdisplayed), continue; end
  LCL = CAND(:,util.X(iindividual,:)==1);  
  lr = fitglm(LCL(li.trnindices,:),li.trnlabels,'linear','Distribution','binomial','BinomialSize',1,Link='logit');
  radscores = LCL*lr.Coefficients.Estimate(2:end) + lr.Coefficients.Estimate(1);
  [radscores,param.radiomics.radscoremus,param.radiomics.radscoresigmas] = stdfeatures(radscores,li.trnindices);
  
  NEW = [radscores CLN];
  [newsubset,B2,fitinfo2] = fs4new(NEW(li.trnindices,:),li.trnlabels,li.cvcvp,param.std.alpha,param.std.clnoption);  
  TRN = NEW(li.trnindices,newsubset);
  [nom,cvscores] = lrm(TRN,li.trnlabels,param.std.hasintercept,li.cvcvp,true);
  [mcvauc,aucs] = meancvaucs(li.cvcvp,li.trnlabels,cvscores);
  nomcoefs = nom.Coefficients.Estimate; if param.std.hasintercept, nomcoefs=nomcoefs(2:end); end
  trnscores = cpredict(nom,TRN,param.clnames);

  % evalua nomogram
  trnresult.nom = perfresult(li.trnlabels,trnscores,param.performance.version);
  valresult.nom = perfresult(li.trnlabels,cvscores,param.performance.version);
  disp(trnresult.nom);
  disp(valresult.nom);
end