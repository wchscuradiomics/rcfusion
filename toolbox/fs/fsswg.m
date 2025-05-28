function subset = fsswg(TRN,trnlabels,criterion,nsteps,penter)
%subset=fsswg(TRN,labels,criterion,nsteps,penter) selects features with stepwiseglm.

if nargin == 2
  criterion = 'Deviance';
  penter = 0.05;
  nsteps = 200;
elseif nargin == 3
  nsteps = 200;
end

if (strcmpi(criterion,'deviance') || strcmpi(criterion,'sse')) && ...
  (nargin==4 || isempty(penter))
  penter = 0.05;
elseif (strcmpi(criterion,'aic') || strcmpi(criterion,'bic') || ...
    strcmpi(criterion,'AdjRsquared')) && (nargin==4 || isempty(penter))
  penter = 0;
elseif strcmpi(criterion,'Rsquared') && (nargin==4 || isempty(penter))
  penter = 0.1;
end

% 'Deviance' 'SSE' 'AIC' 'BIC' 'Rsquared' 'AdjRsquared'
stepmdl = stepwiseglm(TRN,trnlabels,'constant', ...
  'Distribution','binomial', ...
  'Link','logit',...
  'NSteps',nsteps, ...
  'Verbose',0, ...
  'Criterion',criterion, ...
  'pEnter',penter);
subset = mdl2subset(stepmdl); 