function [nom,cvscores] = lrm(D,labels,intercept,cvcvp,isrtnmodel)
%[nom,cvscores]=lrm(D,labels,intercept,cvcvp) fits a linear regression model. If argument cvcvp exists, then cross-validation is
%performed and returned argument cvscores represents cross-validation scores.

cvscores = zeros(height(D),1);
if nargin >= 4
  for k=1:cvcvp.NumTestSets
    mdl = fitglm(D(cvcvp.training(k),:),labels(cvcvp.training(k)),'linear','Distribution','binomial','BinomialSize',1, ...
      Link='logit',Intercept=intercept);
    cvscores(cvcvp.test(k)) = predict(mdl,D(cvcvp.test(k),:));
  end
end

if ~exist("isrtnmodel","var"), isrtnmodel=false; end
if nargin == 3 || isrtnmodel
  nom = fitglm(D,labels,'linear','Distribution','binomial','BinomialSize',1,Link='logit',Intercept=intercept);
else
  nom = [];
end