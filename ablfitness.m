function fitnesses = ablfitness(X,RAD,CLN,labels,cvcvp,util) % used for GLSP (ablation experiment)
%fitnesses=ablfitness(X,RAD,CLN,labels,cvcvp,util) calculates fitness for population X. 
%
% Order of processes: creation, fitness, selection, crossover, muation, check for termination, fitness, ....
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
% szlocal = util.param.abl.szlocal; % size of the local space
idfmax = util.param.abl.idfmax;
generation = util.generation + 1;
clnoption = util.param.abl.clnoption;
alpha = util.param.abl.alpha;
hasintercept = util.param.abl.hasintercept;
popsize = height(X);
fitnesses = ones(popsize,1);
contributionss = cell(height(X),1);
lclsubsets = cell(height(X),1);
[nsmps,nradvars] = size(RAD); 
% nclnvars = width(CLN); 

parfor i=1:height(X) % variables for parallel computing -> suffix '4par'
  warning('off');
  % initialize
  individual = X(i,:); % the i-th individual
  contributions = zeros(1,length(individual)); contributionss{i} = contributions; % default contribution is zero
  rad1indices = individual==1; 
  RAD4PAR = RAD;
  LCL = RAD4PAR(:,rad1indices); 
  % nlclvars = width(LCL);

  % selecting/distilling features from radiomic features
  [lclsubset,B,fitinfo,bstlindex] = fs4lcl(LCL,labels,cvcvp,alpha,fix(idfmax*nsmps/10));
  coefs = B(:,bstlindex); intercept = fitinfo.Intercept(bstlindex);
  radscores = LCL*coefs + intercept;
  [radscores,mus,sigmas] = stdfeatures(radscores,(1:height(radscores))'); 

  % selecting input variables to construct a nomogram and calculating contribution and fitness
  NEW = [radscores CLN];
  [newsubset,B2,fitinfo2] = fs4new(NEW,labels,cvcvp,alpha,clnoption);
  if newsubset(1)~=1 || length(lclsubset) < 4, continue; end  
  [nom,cvscores] = lrm(NEW(:,newsubset),labels,hasintercept,cvcvp,true);
  % [~,~,~,cvauc] = perfcurve(labels,cvscores,1);
  [mcvauc,aucs] = meancvaucs(cvcvp,labels,cvscores);  

  % vc = std(aucs)/mcvauc; if vc > 0.15, continue; end  
  lclsubsets{i} = lclsubset;
  fitnesses(i) = -mcvauc; 
end

util.contributionss = contributionss;
util.fitnesses = fitnesses;
util.lclsubsets = lclsubsets;
util.generation = generation;