function [trnresult,valresult,lclsubset,newsubset,param] = lassonomogram(CAND,CLN,li,param,prelexptparam)

% select radiomic features based on lasso and build a radscore
% rng(param.seed);
[lclsubset,~,B,fitinfo] = fslasso(CAND(li.trnindices,:),li.trnlabels, li.cvcvp,'none|nan', ...
    true,prelexptparam.criteria,prelexptparam.alpha,false,prelexptparam.dfmax);
bstlindex = fitinfo.bestlambdaindex;
coefs = B(:,bstlindex); intercept = fitinfo.Intercept(bstlindex);
radcvscores = fitinfo.YFIT(:,bstlindex);
radscores = CAND*coefs + intercept;
[radmcvauc,radkcvaucs] = meancvaucs(li.cvcvp,li.trnlabels,radcvscores);
[radscores,param.radiomics.radscoremus,param.radiomics.radscoresigmas] = stdfeatures(radscores(:,1),li.trnindices);
% LCL = CAND(:,lclsubset);

% results for radiomics
radtrnscores = radscores(li.trnindices,:);
trnresult.rad = perfresult(li.trnlabels,radtrnscores,param.performance.version);
valresult.rad = perfresult(li.trnlabels,radcvscores,param.performance.version);
valresult.rad.mcvauc = radmcvauc; valresult.rad.stdmcvauc = std(radkcvaucs,0);

% build a nomogram based on NEW
NEW = [radscores CLN];
param.nom.fsmethod = prelexptparam.nomfsmethod;
if strcmpi(prelexptparam.nomfsmethod,'lasso')
  [newsubset,B2,fitinfo2] = fs4new(NEW(li.trnindices,:),li.trnlabels,li.cvcvp,param.lasso.alpha,param.lasso.clnoption);
else
  param.nom.pthreshold = str2double(parstr(prelexptparam.nomfsmethod,2,'|'));  
  [param.nom.pvalues,~,newsubset] = univaranalysis(NEW(li.trnindices,:),li.trnlabels,'regression',param.nom.pthreshold);
end
TRN = NEW(li.trnindices,newsubset);
[nom,nomcvscores] = lrm(TRN,li.trnlabels,param.lasso.hasintercept,li.cvcvp,true);
[nommcvauc,nomkcvaucs] = meancvaucs(li.cvcvp,li.trnlabels,nomcvscores);
nomcoefs = nom.Coefficients.Estimate; if param.lasso.hasintercept, nomcoefs=nomcoefs(2:end); end
radratio = abs(nomcoefs(1))/sum(abs(nomcoefs));
nomtrnscores = cpredict(nom,TRN,param.clnames);

% results for nomogram
trnresult.nom = perfresult(li.trnlabels,nomtrnscores,param.performance.version);
valresult.nom = perfresult(li.trnlabels,nomcvscores,param.performance.version);
valresult.nom.mcvauc = nommcvauc; valresult.nom.stdmcvauc = std(nomkcvaucs);