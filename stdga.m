function [x,util] = stdga(RAD,CLN,labels,cvcvp,param)

nvars = width(RAD);
popsize = 200; % population size
elitecount = ceil(0.05*popsize);
xoverfraction = 0.85; % (popsize - elitecount)/(popsize-elitecount);
clear util;
util = OurUtilHandle;
util.param = param;
util.generation = 0;

options = optimoptions('ga', ... % gamultiobj ga
  'PopulationType', 'bitstring', ...
  'PopulationSize', popsize, ...
  'UseVectorized', true, ...
  'MaxGenerations', 126, ...
  'PlotFcn', 'gaplotbestf', ... % gaplotpareto gaplotbestf gaplotspread {[]} gaplotrankhist  
  'FunctionTolerance', 1e-6, ...  
  ... % The number of selected radiomics features must be limited to reduce overfitting. % Worse performance with ourcreation
  'CreationFcn', @(nvars,f,options)ourcreation(options.PopulationSize,RAD,CLN,labels,cvcvp,util), ... 
  'EliteCount', elitecount, ...
  'CrossoverFraction', xoverfraction, ...  
  'display','off' ...
  ); 
  
x=ga(@(X)stdfitness(X,RAD,CLN,labels,cvcvp,util),nvars,[],[],[],[],[],[],[],options);