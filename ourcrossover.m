function xoverkids=ourcrossover(parents,POP,util)
%xoverkids=ourcrossover
%
% POP: a popsize-by-nvars matrix representing this (the current) population.
%
% % how many crossover offspring will there be from each source?
% nEliteKids = options.EliteCount;
% nXoverKids = round(options.CrossoverFraction * (size(thisPopulation,1) - nEliteKids));
% nMutateKids = size(thisPopulation,1) - nEliteKids - nXoverKids;
% % how many parents will we need to complete the population?
% nParents = 2 * nXoverKids + nMutateKids;
%
% Randomly pair up two individuals from the parents to produce offspring, requiring that the two individuals cannot be the same
% (the presence of repeated indexes in the parents indicates that they can team up with multiple other individuals to produce
% offspring);

if ~isequal(util.X,POP), error('No consistent Population.'); end
nvars = width(POP); % szlocal = util.param.our.szlocal;
xoverkids = zeros(length(parents)/2,nvars);

% step 1/2: select parents by pairwise
sn = parents(randperm(length(parents)));  
while ~chkrndparents(sn), sn = parents(randperm(length(parents))); end

% setp 2/2: a new individual takes half of their genes from each parent
k=1;
for i=1:2:length(sn)    
  xoverkids(k,:)=crossover(POP(sn(i),:),POP(sn(i+1),:),sn(i),sn(i+1),util); 
  k = k+1;
end
end % end for function ourcrossover

%% check sn(i) and sn(i+1) for a pairwise crossover
function r=chkrndparents(sn)
A = zeros(length(sn)/2,2);
k = 0;
r = true; 
for i=1:2:length(sn)
  k = k + 1;
  A(k,:) = sort([sn(i) sn(i+1)]);
  if sn(i)==sn(i+1)
    r=false;
  end
end
if ~r, return; end
B = unique(A,'rows');
if height(A)~= height(B), r = false; end
end

%% create a new individual
function offspring=crossover(a,b,ia,ib,util)

nvars = length(a); % szlocal = util.param.our.szlocal;

% check a, b, ia, ib, util.X
if ~(isequal(a,util.X(ia,:)) && isequal(b,util.X(ib,:))), error('a and b are not consistent with util.X(ia and ib,:)'); end

conmax = max([util.contributionss{ia} util.contributionss{ib}]);

offspring = zeros(1,nvars);
for j=1:nvars
  if a(j) == 0 && b(j) == 0
  elseif a(j)==1 && b(j) == 1
    offspring(j) = 1;
  elseif a(j) == 0 && b(j) == 1
    con = util.contributionss{ib}(j);
    if rand <= con/conmax, offspring(j) = 1; end
  elseif a(j) == 1 && b(j) == 0
    con = util.contributionss{ia}(j);
    if rand <= con/conmax, offspring(j) = 1; end
  end
end

while sum(offspring) == 0 % use the default crossover method
  selector = rand(1,nvars) <= 0.5;
  aindices = find(selector); bindices = find(~selector);
  offspring(aindices) = a(aindices);
  offspring(bindices) = a(bindices);
end

end % end function crossover 