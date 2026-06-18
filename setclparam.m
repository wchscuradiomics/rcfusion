function param = setclparam(mls,opt)

if nargin==1, opt.suffix=''; opt.GE=[]; opt.gelabels=[]; opt.ES=[]; opt.eslabels=[]; end

if ~isfield(opt,'GE'), opt.GE=[]; end
if ~isfield(opt,'gelabels'), opt.gelabels=[]; end
if ~isfield(opt,'ES'), opt.ES=[]; end
if ~isfield(opt,'eslabels'), opt.eslabels=[]; end

param.comoption = gsdlftcomoption(); 
param.comoption.gerate = 0.1;
param.comoption.optobj = 'max_mcvauc';
param.comoption.gestrategy = 'geoff';
param.comoption.esstrategy = 'noes';
param.comoption.cloptimizer = 'cstmgridsearch'; % bi

for iml = 1:length(mls) % param.(ml) is an option for ml
  ml = lower(mls{iml});
  param.(ml) = param.comoption; % common option for all ML approaches
  if strcmpi(param.(ml).gestrategy,'geon'), param.(ml).GE = opt.GE; param.(ml).gelabels = opt.gelabels; end

  switch ml
    case 'lr'
      param.lr.useintercept = true;
    case 'svm'
      param.svm.kernel = 'all';  
    case 'nb'
      param.nb.vartypes = zeros(1,width(opt.varcount));
    case 'lgb'
      param.lgb.varpath = ['lgb_' opt.suffix '.pkl'];
      switch param.lgb.esstrategy
        case {'noes','esincv'} % (defult) use n_estimators as a hyperparameter         
          param.lgb.ES = [];
          param.lgb.eslabels = [];
        case 'eson' % use an additional validation set for early stopping
          param.lgb.ES = opt.ES;
          param.lgb.eslabels = opt.eslabels;
      end
    case 'xgb'
      param.xgb.varpath = ['xgb_' opt.suffix '.pkl'];
      switch param.xgb.esstrategy
        case {'noes','esincv'} % (defult) use n_estimators as a hyperparameter 
          param.xgb.ES = [];
          param.xgb.eslabels = [];
        case 'eson' % use an additional validation set for early stopping
          param.xgb.ES = opt.ES;
          param.xgb.eslabels = opt.eslabels;
      end
    case 'ann'
      switch param.ann.esstrategy
        case {'noes','esincv'} % (defult) use n_estimators as a hyperparameter
          param.ann.ES = [];
          param.ann.eslabels = [];
        case 'eson' % use an additional validation set for early stopping
          param.ann.ES = opt.ES;
          param.ann.eslabels = opt.eslabels;
      end
    case {'adaboostm1', 'logitboost', 'gentleboost', 'robustboost', 'lpboost', 'totalboost', 'rusboost'}
      param.(ml).ml = ml;
    otherwise
  end

end % end for each ML approach