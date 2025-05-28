function [mcvauc,aucs] = meancvaucs(cvcvp,labels,cvscores)
%mcvauc=meancvaucs(cvcvp,cvscores) calculates a mean CV (cross-validation)
%AUC for k-fold scores.
%
% cvcvp: cross-validation_cvpartition, a CV partition object representing
% k-fold cross-validation.
%
% labels: a n-by-1 vector of double type representing labels of n samples,
% where 1 representing a positive class.
% 
% cvscores: cross-validation_scores, a n-by-1 or n-by-m matrix representing
% scores of n samples, of which the j-th column corresponing to the scores
% of the j-th class and the first column is used for the positive class of
% function perfcurve.

if ~isequal(unique(labels),[0;1]), error('labels ~= 0|1 classes.'); end

cvscores = cvscores(:,1);
aucs = NaN(1, cvcvp.NumTestSets);
for j=1:cvcvp.NumTestSets
  [~,~,~,aucs(j)] = perfcurve(labels(test(cvcvp,j)),...
    cvscores(test(cvcvp,j)),1);
end
mcvauc = mean(aucs,'all','includenan');