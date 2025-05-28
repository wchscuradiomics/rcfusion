function SC = cpredict(model,D,clnames)
%SC=cpredict(model,D,clnames) predicts scores for D. SC is a n-by-2 matrix, where n is the number of samples and the first column
%represents the positive class.

if nargin == 2, clnames = [1 0]; end

if isa(model,'GeneralizedLinearModel')
  if ~strcmpi(model.Distribution.Name, 'Binomial') || ~strcmpi(model.Link.Name,'logit')
    error('Just the distribution of binomial with size 1 and the link of ligit are supported.');
  end
  % disp('Please insure the distribution is binomial with size 1 and the link is logit.');
  slr = predict(model,D);
  SC = [slr 1-slr];
elseif isa(model,'network')
  SC = model(D')'; % model.pnn(DS')' or sim(model.pnn,DS')'
elseif isa(model,'py.lightgbm.sklearn.LGBMClassifier')
  SC = lgbpredict(model,D,clnames);
elseif isa(model,'py.xgboost.sklearn.XGBClassifier')
  SC = xgbpredict(model,D,clnames);
elseif isfield(model,'name') &&  strcmpi(model.name,'lasso')
  if isfield(model,'coefs'), coefs = model.coefs; else, coefs = model.coef; end
  if width(D) ~= length(coefs), coefs(abs(coefs)<=1e-10) = []; end  
  yhats = D*coefs + model.intercept;  
  SC = [yhats -yhats];
else
  [~,SC] = predict(model,D);
end