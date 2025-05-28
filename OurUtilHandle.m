classdef OurUtilHandle < handle % used for GLoRC, GLSP, FCM, RPP, and GA
  properties
    % Note: 
    % 1. popsize i.e., population size
    % 2. nvars represents the number of radiomic features

    X % population
    INIT % initial population
    contributionss % a popsize-by-1 cell representing contributions of popsize individuals
    param % settings for selecting and evaluating features
    fitnesses % a popsize-by-1 array representing fitnesses of popsize individuals
    lclsubsets % local subset, i.e., indices of distilled features
    generation % representing the current generation
  end % end for properties

  methods
    
  end % end for methods
end % end for classdef