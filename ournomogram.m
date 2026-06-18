function [x,util] = ournomogram(CAND,CLN,li,param)

% subset based on ga
rng(param.seed);
[x,util] = ourga(CAND(li.trnindices,:),CLN(li.trnindices,:),li.trnlabels,li.cvcvp,param); % only training samples are used

% use x to construct a nomogram (only on training samples) (an independent test set in CAND is used to test)
for iindividual = 1:height(util.X)
  if ~isequal(util.X(iindividual,:),x),continue; end

  % use the subset of radiomics features to compute a rad-score (only on training samples)
  LCL = CAND(:,util.X(iindividual,:)==1);  
  [lclsubset,B,fitinfo,bstlindex] = fs4lcl(LCL(li.trnindices,:),li.trnlabels,li.cvcvp,param.our.alpha); % using training samples
  coefs = B(:,bstlindex); intercept = fitinfo.Intercept(bstlindex);
  radscores = LCL*coefs + intercept;
  [radscores,param.radiomics.radscoremus,param.radiomics.radscoresigmas] = stdfeatures(radscores,li.trnindices);
  
  % construct a nomogram (only on training samples)
  NEW = [radscores CLN];
  [newsubset,B2,fitinfo2] = fs4new(NEW(li.trnindices,:),li.trnlabels,li.cvcvp,param.our.alpha,param.our.clnoption); % using training samples 
  TRN = NEW(li.trnindices,newsubset);
  [nom,cvscores] = lrm(TRN,li.trnlabels,param.our.hasintercept,li.cvcvp,true);  % using training samples
  [mcvauc,aucs] = meancvaucs(li.cvcvp,li.trnlabels,cvscores);
  nomcoefs = nom.Coefficients.Estimate; if param.our.hasintercept, nomcoefs=nomcoefs(2:end); end  
  
  % validate and test the nomogram
  trnscores = cpredict(nom,TRN,param.clnames);
  TST = NEW(li.tstindices,newsubset);  
  tstscores = cpredict(nom,TST,param.clnames);
  trnresult.nom = perfresult(li.trnlabels,trnscores,param.performance.version);
  valresult.nom = perfresult(li.trnlabels,cvscores,param.performance.version);
  tstresult.nom = perfresult(li.tstlabels,tstscores,param.performance.version);

  % display result
  s = [...
    num2str(util.generation) 9 num2str(trnresult.nom.auc(1)) 9 num2str(valresult.nom.auc(1)) 9 num2str(tstresult.nom.auc(1)) 9 ...
    num2str(abs(nomcoefs(1))/sum(abs(nomcoefs))) 9 num2str(1-length(util.lclsubsets{iindividual})/width(CAND)) 9 ...
    num2str(util.fitnesses(iindividual)) 9 num2str(sum(util.X(iindividual,:))) 9 ...
    num2str(length(util.lclsubsets{iindividual})) 9 num2str(length(newsubset)) ...
  ];
  disp(['x=' num2str(iindividual) 9 s]); 
  break;  
end