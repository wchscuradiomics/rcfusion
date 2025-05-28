function fitnesses = ourfitness(X,RAD,CLN,labels,cvcvp,util)
%fitnesses=ourfitness(X,RAD,CLN,labels,cvcvp,util) calculates fitness for population X.
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
szlocal = util.param.our.szlocal;
generation = util.generation + 1;
radconrithreshold = util.param.our.radconrithreshold;
clnoption = util.param.our.clnoption;
alpha = util.param.our.alpha;
hasintercept = util.param.our.hasintercept;
fitcoefs = util.param.our.fitcoefs;
popsize = height(X);
fitnesses = ones(popsize,1);
contributionss = cell(height(X),1);
lclsubsets = cell(height(X),1);
[nsmps,nradvars] = size(RAD); nclnvars = width(CLN); 

parfor i=1:height(X) % variables for parallel computing -> suffix '4par'
  warning('off');
  % initialize
  individual = X(i,:); % the i-th individual
  contributions = zeros(1,length(individual)); contributionss{i} = contributions; % default contribution is zero
  rad1indices = find(individual==1); 
  LCL = RAD(:,rad1indices); nlclvars = width(LCL);
  abc = fitcoefs;

  % LASSO for selecting features from radiomics and clinical features
  [lclsubset,B,fitinfo,bstlindex] = fs4lcl(LCL,labels,cvcvp,alpha);
  coefs = B(:,bstlindex); intercept = fitinfo.Intercept(bstlindex);
  radscores = LCL*coefs + intercept;
  [radscores,mus,sigmas] = stdfeatures(radscores,(1:height(radscores))');  
  NEW = [radscores CLN];
  [newsubset,B2,fitinfo2] = fs4new(NEW,labels,cvcvp,alpha,clnoption);
 
  % calculate contribution and fitness
  if newsubset(1) ~= 1 || length(lclsubset) < 4, continue; end  
  [nom,cvscores] = lrm(NEW(:,newsubset),labels,hasintercept,cvcvp,true);
  % [~,~,~,cvauc] = perfcurve(labels,cvscores,1);
  [mcvauc,aucs] = meancvaucs(cvcvp,labels,cvscores);  

  % vc = std(aucs)/mcvauc; if vc > 0.15, continue; end
  nomcoefs = nom.Coefficients.Estimate; if hasintercept, nomcoefs=nomcoefs(2:end); end
  contributions(rad1indices(lclsubset)) = abs(coefs(lclsubset))/sum(abs(coefs(lclsubset)));
  contributionss{i} = contributions*mcvauc*abs(nomcoefs(1))/sum(abs(nomcoefs));
  radconratio = abs(nomcoefs(1))/sum(abs(nomcoefs));  
  lclsellenratio = 1-length(lclsubset)/nradvars; % fix(nsmps/10);
  if radconratio > radconrithreshold, radconratio=radconrithreshold-radconratio+radconrithreshold; end
  lclsubsets{i} = lclsubset;
  fitnesses(i) = -( mcvauc*abc(1)+radconratio*abc(2) + abc(3)*lclsellenratio );  
end

util.contributionss = contributionss;
util.fitnesses = fitnesses;
util.lclsubsets = lclsubsets;
util.generation = generation;