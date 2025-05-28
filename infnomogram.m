function [trnresult,valresult,lclsubset,newsubset,param] = infnomogram(CAND,CLN,li,param,prelexptparam)

% select radiomic features based on il and build a radscore
% rng(param.seed);
if strcmpi(prelexptparam.way, 'U')
  [radfeatsindices, ~, lclsubset] = InfFS_U( CAND(li.trnindices,:), li.trnlabels, prelexptparam.alpha1, 0 );
elseif strcmpi(prelexptparam.way, 'S')
  [radfeatsindices, ~, lclsubset] = InfFS_S(CAND(li.trnindices,:), li.trnlabels , ...
    [prelexptparam.alpha1 prelexptparam.alpha2 prelexptparam.alpha3]);
else
  error(['Not supported way: ' prelexptparam.way]);
end
if length(lclsubset)<4, lclsubset=radfeatsindices(1:4); end

LCL = CAND(:,lclsubset);
[lr,radcvscores] = lrm(LCL(li.trnindices,:),li.trnlabels,true,li.cvcvp,true);
coefs = lr.Coefficients.Estimate(2:end); intercept = lr.Coefficients.Estimate(1);
radscores = LCL*coefs + intercept;
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
  [newsubset,B2,fitinfo2] = fs4new(NEW(li.trnindices,:),li.trnlabels,li.cvcvp,param.inf.alpha,param.inf.clnoption);
elseif contains(prelexptparam.nomfsmethod,'univaranalysis','IgnoreCase',true)
  param.nom.pthreshold = str2double(parstr(prelexptparam.nomfsmethod,2,'|'));  
  [param.nom.pvalues,~,newsubset] = univaranalysis(NEW(li.trnindices,:),li.trnlabels,'regression',param.nom.pthreshold);
else % inf
  param.nom.way = prelexptparam.way;
  param.nom.alpha1 = prelexptparam.alpha1; param.nom.alpha2 = prelexptparam.alpha2; param.nom.alpha3 = prelexptparam.alpha3;
  if strcmpi(prelexptparam.way, 'U')
    [nomfeatsindices, ~, newsubset] = InfFS_U(NEW(li.trnindices,:),li.trnlabels,prelexptparam.alpha1,0);
  else
    [nomfeatsindices, ~, newsubset] = InfFS_S(NEW(li.trnindices,:),li.trnlabels, ...
      [prelexptparam.alpha1 prelexptparam.alpha2 prelexptparam.alpha3]);
  end  
  % To improve performance, at least 3 clinical variables are included.
  if length(newsubset) < 4, newsubset = nomfeatsindices(1:4); end 
end
TRN = NEW(li.trnindices,newsubset);
[nom,nomcvscores] = lrm(TRN,li.trnlabels,param.inf.hasintercept,li.cvcvp,true);
[nommcvauc,nomkcvaucs] = meancvaucs(li.cvcvp,li.trnlabels,nomcvscores);
nomcoefs = nom.Coefficients.Estimate; if param.inf.hasintercept, nomcoefs=nomcoefs(2:end); end
radratio = abs(nomcoefs(1))/sum(abs(nomcoefs));
nomtrnscores = cpredict(nom,TRN,param.clnames);

% results for nomogram
trnresult.nom = perfresult(li.trnlabels,nomtrnscores,param.performance.version);
valresult.nom = perfresult(li.trnlabels,nomcvscores,param.performance.version);
valresult.nom.mcvauc = nommcvauc; valresult.nom.stdmcvauc = std(nomkcvaucs);