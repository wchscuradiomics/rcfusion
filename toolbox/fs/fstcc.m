function subset = fstcc(TRN,tralabels,option,pthreshold,cthreshold)
%subset = fstcc(TRN,tralabels,option,pthreshold,cthreshold) perform feature selection using two stages feature selection using two
%stages comprised of significance difference tests and Fisher's correlation coefficients.
%
% DISCARDED, Please use fs3multistages instead.
%
% Currently, the significance test of difference used the t-test method, the default threshold p-value is set to 0.05, and the
% correlation threshold is set to 0.8.
%
% A row vector is a sample.
%
% option: 'two-stages' (default) /'two-stages:utest'/'two-stages:ttest' | 'correlation' | 'ttest' | 'utest' .

warning('DISCARDED, Please use fs3multistages instead.');

names = unique(tralabels);
% if length(names) ~=2, error('There must be two types in labels.'); end

if nargin == 2
  option = 'two-stages';
  pthreshold=0.05; cthreshold=0.8;
elseif nargin==3
  pthreshold=0.05; cthreshold=0.8;
end

%% 开始预选，分3种情况：
if contains(option,'two-stages')
  if contains(option,'utest')
    [p,~] = mwu(TRN(tralabels==names(1),:),TRN(tralabels==names(2),:)); % u test
  else
    [~,p] = ttest2(TRN(tralabels==names(1),:),TRN(tralabels==names(2),:));
  end
  presubset = find(p<pthreshold);
  further = fscc(TRN(:,presubset),cthreshold);
  subset = presubset(further);
elseif strcmpi(option,'correlation')  
  subset = fscc(TRN,cthreshold);
elseif strcmpi(option,'ttest')
  [~,p] = ttest2(TRN(tralabels==names(1),:),TRN(tralabels==names(2),:));
  subset = find(p<pthreshold);
elseif strcmpi(option,'utest')
  [p,~] = mwu(TRN(tralabels==names(1),:),TRN(tralabels==names(2),:)); % u test
  subset = find(p<pthreshold);
else
  error('Invalid option. The parameter option can be ''two-stages'' (default) | ''correlation'' | ''ttest'' ');
end