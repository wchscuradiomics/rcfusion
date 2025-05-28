function [trnresult,valresult,lclsubset,newsubset,param] = rfenomogram(CAND,CLN,li,param,prelexptparam)
%%% comment it as a script for final experiment

% select radiomic features based on rfe and build a radscore
% rng(param.seed);
% lr = stepwiseglm(CAND(li.trnindices,:),li.trnlabels,'constant','Distribution','binomial','Upper','linear','Link','logit',...
%    'NSteps',width(CAND),'Verbose',0,'Criterion',prelexptparam.criterion); 
% lclsubset = mdl2subset(lr);
[lclsubset,~] = fsrfe(CAND(li.trnindices,:),li.trnlabels,li.cvcvp,prelexptparam.nlim,'glm',param.seed,false);
if length(lclsubset) < 3, trnresult=[]; valresult=[]; lclsubset=[]; newsubset=[]; return; end
LCL = CAND(:,lclsubset);
[lr,radcvscores] = lrm(LCL(li.trnindices,:),li.trnlabels,true,li.cvcvp,true);
coefs = lr.Coefficients.Estimate(2:end); intercept = lr.Coefficients.Estimate(1);
try
  radscores = LCL*coefs + intercept;
catch ex
  disp(lr);
  disp(size(LCL));
  disp(size(coefs));
  error(ex.message);
end
[radmcvauc,radkcvaucs] = meancvaucs(li.cvcvp,li.trnlabels,radcvscores);
[radscores,param.radiomics.radscoremus,param.radiomics.radscoresigmas] = stdfeatures(radscores(:,1),li.trnindices);

% results for radiomics
radtrnscores = radscores(li.trnindices,:);
trnresult.rad = perfresult(li.trnlabels,radtrnscores,param.performance.version);
valresult.rad = perfresult(li.trnlabels,radcvscores,param.performance.version);
valresult.rad.mcvauc = radmcvauc; valresult.rad.stdmcvauc = std(radkcvaucs,0);

% build a nomogram based on NEW
NEW = [radscores CLN];
param.nom.fsmethod = prelexptparam.nomfsmethod;
if strcmpi(prelexptparam.nomfsmethod,'lasso')
  [newsubset,B2,fitinfo2] = fs4new(NEW(li.trnindices,:),li.trnlabels,li.cvcvp,param.rfe.alpha,param.rfe.clnoption);
elseif strcmpi(prelexptparam.nomfsmethod,'rfe')
  [newsubset,history] = fsrfe(NEW(li.trnindices,:),li.trnlabels,li.cvcvp,param.rfe.nlim,'glm',param.seed,false);
elseif contains(prelexptparam.nomfsmethod,'univaranalysis','IgnoreCase',true)
  param.nom.pthreshold = str2double(parstr(prelexptparam.nomfsmethod,2,'|'));  
  [param.nom.pvalues,~,newsubset] = univaranalysis(NEW(li.trnindices,:),li.trnlabels,'regression',param.nom.pthreshold);
else 
  error(['Not supported method: ' prelexptparam.nomfsmethod]);
end
TRN = NEW(li.trnindices,newsubset);
[nom,nomcvscores] = lrm(TRN,li.trnlabels,param.rfe.hasintercept,li.cvcvp,true);
[nommcvauc,nomkcvaucs] = meancvaucs(li.cvcvp,li.trnlabels,nomcvscores);
nomcoefs = nom.Coefficients.Estimate; if param.rfe.hasintercept, nomcoefs=nomcoefs(2:end); end
radratio = abs(nomcoefs(1))/sum(abs(nomcoefs));
nomtrnscores = cpredict(nom,TRN,param.clnames);

% results for nomogram
trnresult.nom = perfresult(li.trnlabels,nomtrnscores,param.performance.version);
valresult.nom = perfresult(li.trnlabels,nomcvscores,param.performance.version);
valresult.nom.mcvauc = nommcvauc; valresult.nom.stdmcvauc = std(nomkcvaucs);