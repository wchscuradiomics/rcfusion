function fltrresult = fltrprelexptresult(T,ml,lennomthreshold,mode,lclsellenratio)
%fltrresult=fltrprelexptresult(T,ml,mode,lclsellenratio) filters preliminary experiment results.
%
% T: a n-by-m table representing preliminary experiment results.
%
% ml: a character vector representing the building method based on linear machine learning approaches (Logistic Regression,
% regularization Regression (LASSO, Elastic Net, etc.)
%
% mode: 'complete_independence' (default) | 'on_test'. 'complete_independence' represents selectiing a optmized mode/result from T
% based on training results (including cross-validation results on the training set) and it is completely independent of the test.
% 'on_test' represents preliminary experiment results are based on the training set, but it selects a ml-based model based on test
% results.

if nargin==3, mode='complete_independence'; end

mcvaucs = T.([ml 'mcvauc']);
cvaucs = T.([ml 'cvauc']);
cvaccs = T.([ml 'cvacc']);
tstaucs = T.([ml 'tstauc']);
tstaccs = T.([ml 'tstacc']);
trnaucs = T.([ml 'trnauc']);
stdmcvaucs = T.([ml 'stdmcvauc']); 
lenrads = T.('lenradscores');
lennoms = T.('lennoms');

% trnaucs < cvaucs or trnaucs < mcvaucs means a potentially bad fit.
if strcmpi(mode,'complete_independence')
  vldindices = find(trnaucs > cvaucs & trnaucs > mcvaucs & lennoms <= lennomthreshold); %  
  maxobj = max(mcvaucs(vldindices));
  vldindices = vldindices(mcvaucs(vldindices) >= maxobj*0.95); % maxobj - mcvaucs(vldindices) <= 0.01  
  [~, svldindices] = max((1-lclsellenratio)*mcvaucs(vldindices) - lclsellenratio*lenrads(vldindices)/51);
  if isempty(svldindices), fltrresult=[]; return; end
  fltrresult = T(vldindices(svldindices(1)),:);
elseif strcmpi(mode,'on_test')
  vldindices = find(trnaucs > cvaucs & trnaucs > mcvaucs & lennoms <= lennomthreshold ... % bad fits are also removed
    & mcvaucs > tstaucs & cvaucs > tstaucs & cvaccs > tstaccs); 
  maxobj = max(tstaucs(vldindices));
  vldindices = vldindices(tstaucs(vldindices) >= maxobj*0.95 );
  [~, svldindices] = max((1-lclsellenratio)*tstaucs(vldindices) - lclsellenratio*lenrads(vldindices)/51);  
  if isempty(svldindices), fltrresult=[]; return; end
  fltrresult = T(vldindices(svldindices(1)),:);
end


