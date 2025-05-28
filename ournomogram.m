function [x,util] = ournomogram(CAND,CLN,li,param)

% subset based on ga
rng(param.seed);
[x,util] = ourga(CAND(li.trnindices,:),CLN(li.trnindices,:),li.trnlabels,li.cvcvp,param); 
isxdisplayed = false;
for iindividual = 1:height(util.X)
  if (isequal(util.X(iindividual,:),x) && isxdisplayed) || ~isequal(util.X(iindividual,:),x),continue; end
  % if (isequal(util.X(iindividual,:),x) && isxdisplayed), continue; end
  LCL = CAND(:,util.X(iindividual,:)==1);  
  [lclsubset,B,fitinfo,bstlindex] = fs4lcl(LCL(li.trnindices,:),li.trnlabels,li.cvcvp,param.our.alpha);
  coefs = B(:,bstlindex); intercept = fitinfo.Intercept(bstlindex);
  radscores = LCL*coefs + intercept;
  [radscores,param.radiomics.radscoremus,param.radiomics.radscoresigmas] = stdfeatures(radscores,li.trnindices);
  
  NEW = [radscores CLN];
  [newsubset,B2,fitinfo2] = fs4new(NEW(li.trnindices,:),li.trnlabels,li.cvcvp,param.our.alpha,param.our.clnoption);  
  TRN = NEW(li.trnindices,newsubset);
  [nom,cvscores] = lrm(TRN,li.trnlabels,param.our.hasintercept,li.cvcvp,true);
  [mcvauc,aucs] = meancvaucs(li.cvcvp,li.trnlabels,cvscores);
  nomcoefs = nom.Coefficients.Estimate; if param.our.hasintercept, nomcoefs=nomcoefs(2:end); end
  trnscores = cpredict(nom,TRN,param.clnames);

  % evalua nomogram
  trnresult.nom = perfresult(li.trnlabels,trnscores,param.performance.version);
  valresult.nom = perfresult(li.trnlabels,cvscores,param.performance.version);
  disp(trnresult.nom);
  disp(valresult.nom);
end