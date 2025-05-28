function fitnesses = stdfitness(X,RAD,CLN,labels,cvcvp,util)
%fitnesses=stdfitness(X,RAD,CLN,labels,cvcvp,util) calculates fitness for population X.
%
% Order of processes: creation, fitness, selection, crossover, muation, judge for termination, fitness, ....
%
% popsize i.e., population size;
%
% When 'UseVectorized' is true, then X is a popsize-by-nvars matrix ; else X is a 1-by-nvars vector.
%
% util: OurUtilHandle object
%
% RAD: a n-by-m matrix representing candidate radiomics features and imaging parameters, where n is the number of samples, m is
% the number of candidate features and the (m-1)- and m-th features are X-ray Tube Current and KVP values.

util.X = X; % current population X
klcllen = util.param.std.klcllen;
generation = util.generation + 1;
clnoption = util.param.std.clnoption;
alpha = util.param.std.alpha;
hasintercept = util.param.std.hasintercept;
popsize = height(X);
fitnesses = ones(popsize,1);
[nsmps,nradvars] = size(RAD); nclnvars = width(CLN); 

parfor i=1:height(X) % variables for parallel computing -> suffix '4par'
  warning('off');
  % initialize
  individual = X(i,:); % the i-th individual
  contributions = zeros(1,length(individual)); contributionss{i} = contributions; % default contribution is zero
  rad1indices = find(individual==1); 
  LCL = RAD(:,rad1indices); nlclvars = width(LCL);

  % build radscore from radiomics and clinical features
  lr = fitglm(LCL,labels,'linear','Distribution','binomial','BinomialSize',1,Link='logit');
  radscores = LCL*lr.Coefficients.Estimate(2:end) + lr.Coefficients.Estimate(1);
  [radscores,mus,sigmas] = stdfeatures(radscores,(1:height(radscores))');  
  NEW = [radscores CLN];
  [newsubset,B2,fitinfo2] = fs4new(NEW,labels,cvcvp,alpha,clnoption);
 
  % calculate contribution and fitness
  if newsubset(1)~=1, continue; end  
  [nom,cvscores] = lrm(NEW(:,newsubset),labels,hasintercept,cvcvp,true);
  % [~,~,~,cvauc] = perfcurve(labels,cvscores,1);
  [mcvauc,aucs] = meancvaucs(cvcvp,labels,cvscores);  

  % vc = std(aucs)/mcvauc; if vc > 0.15, continue; end
  fitnesses(i) = -(1-klcllen)*mcvauc - klcllen*(1-nlclvars/nradvars);  
  % fitnesses(i) = -mcvauc; 
end

util.fitnesses = fitnesses;
util.generation = generation;