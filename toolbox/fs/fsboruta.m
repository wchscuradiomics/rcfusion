function [subset,sorting]=fsboruta(TRN,trnlabels,alpha,maxdepth,seed)

if nargin==2
  alpha = 0.05; maxdepth = 5; seed = 1;
elseif nargin==3
  maxdepth = 5; seed = 1;
elseif nargin==4
  seed = 1;
end

nfeatures = width(TRN);

% mod = py.importlib.import_module('borutautil');
% py.importlib.reload(mod);

TRN = py.numpy.array(TRN);
trnlabels = py.numpy.array(trnlabels);
pr = py.borutautil.fsboruta(TRN,trnlabels,alpha,maxdepth,seed);
subset = find(logical(pr{1}));
sorting = double(pr{2});
if length(subset) < 5
  imax = max(sorting);
  subset = false(1,nfeatures);
  for i=1:imax
    subset(sorting==i) = true;
    if sum(subset == true) >= 5, break; end
  end
  subset = find(subset);
end
