function [nextScore,nextPopulation,state] = stepGA(thisScore,thisPopulation,options,state,GenomeLength,FitnessFcn)
%STEPGA Moves the genetic algorithm forward by one generation
%   This function is private to GA.

%   Copyright 2003-2021 The MathWorks, Inc.

% how many crossover offspring will there be from each source?
nEliteKids = options.EliteCount;
nXoverKids = round(options.CrossoverFraction * (size(thisPopulation,1) - nEliteKids));
nMutateKids = size(thisPopulation,1) - nEliteKids - nXoverKids;
% how many parents will we need to complete the population?
nParents = 2 * nXoverKids + nMutateKids;

% decide who will contribute to the next generation

% fitness scaling
state.Expectation = feval(options.FitnessScalingFcn,thisScore,nParents,options.FitnessScalingFcnArgs{:});

% selection. parents are indices into thispopulation
parents = feval(options.SelectionFcn,state.Expectation,nParents,options,options.SelectionFcnArgs{:});

% shuffle to prevent locality effects. It is not the responsibility
% if the selection function to return parents in a "good" order so
% we make sure there is a random order here.
parents = parents(randperm(length(parents)));

[~,k] = sort(thisScore);

% Everyones parents are stored here for genealogy display
state.Selection = [k(1:options.EliteCount);parents'];

% here we make all of the members of the next generation
eliteKids  = thisPopulation(k(1:options.EliteCount),:);
xoverKids  = feval(options.CrossoverFcn, parents(1:(2 * nXoverKids)), ... 
    options,GenomeLength,FitnessFcn,thisScore,thisPopulation, ... 
    options.CrossoverFcnArgs{:});
mutateKids = feval(options.MutationFcn,  parents((1 + 2 * nXoverKids):end), ... 
    options,GenomeLength,FitnessFcn,state,thisScore,thisPopulation, ... 
    options.MutationFcnArgs{:});

% group them into the next generation
state.EvalElites = true;
if state.EvalElites % force eval elites
    nextPopulation = [ eliteKids ; xoverKids ; mutateKids ];
else    
    nextPopulation = [ xoverKids ; mutateKids ]; 
end

% Score the population (unique individuals). Using a double cache is ideal
% but using unique for now which is available in MATLAB
if state.HaveDuplicates
    % Uses 1e-12 as absolute tolerance to detect dups
    DS = max(abs(nextPopulation),[],1);
    tol = 1e-12;
   [pop,IA,IC] = uniquetol(nextPopulation,tol, 'DataScale',DS,'ByRows',true);   
else
    pop = nextPopulation;
end

% We want to add the vectorizer if fitness function is NOT vectorized
if strcmpi(options.Vectorized, 'off')
    nextScore = fcnvectorizer(pop,FitnessFcn,1,options.SerialUserFcn,options.ProblemdefOptions.FunOnWorkers);
else
    nextScore = FitnessFcn(pop);
end

% Make sure score is a column vector
nextScore = nextScore(:);
state.FunEval = state.FunEval + size(nextScore,1);

% Verify if we have duplicates in the popultaion and turn off the flag for
% dup evaluation in future.
if state.HaveDuplicates
    nextScore = nextScore(IC,:);
    % Check if we still have dups
    state.HaveDuplicates = numel(IA) ~= numel(IC);
end

if ~state.EvalElites
    % Add elites in the population
    nextPopulation = [ eliteKids ; nextPopulation ];
    nextScore = [thisScore(k(1:options.EliteCount)); nextScore];    
end

% Verify one time if elites needs to be evaluated again.
if state.Generation == 1 && state.EvalElites
    
    delta = thisScore(k(1:options.EliteCount)) - nextScore(1:options.EliteCount);
    
    if ~isempty(delta) && ...       
       max(abs(delta)) <= eps*abs(max(thisScore(k(1:options.EliteCount))))
       % strcmpi(options.Vectorized, 'off') && ...
       state.EvalElites = false;
    end
end

% Update state struct if required
if ~isempty(options.StateUpdateFcn)
    state = options.StateUpdateFcn(state,options);
end