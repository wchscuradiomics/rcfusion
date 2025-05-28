function children = ourmutation(parents,POP,util)

nvars = width(POP);
children = zeros(length(parents),nvars);
for i=1:length(parents)
  children(i,:) = mutation(POP(parents(i),:),parents(i),util);
end

end

function newindividual = mutation(individual,ii,util)

szlocal = util.param.our.szlocal;
lclsubset = util.lclsubsets{ii};

if ~isequal(individual,util.X(ii,:)), error('Inconsistent individual for argument individual and util.X(ia,:)'); end
nvars = width(individual);
contributions = util.contributionss{ii};
[scontributions,sconindices] = sort(contributions,'descend');
vertex = find(scontributions<=1e-6,1,'first') - 4;
grads = abs(diff(scontributions)); grads = normalize(grads(1:vertex),'range');
convexes = abs(diff2(scontributions)); convexes(1)=[]; convexes = normalize(convexes(1:vertex),'range');
[~,pos] = max(grads-convexes);

% mutate
fracindices = 1:nvars; fracindices(sconindices(1:pos)) = []; % not selected for mutation
newindividual = individual;
% newindividual(individual==1 & contributions==0) = 0; % mutate 1 to 0
newindividual( fracindices(individual(fracindices) == 1) ) = 0;
fracindices(individual(fracindices) == 1) = [];
p = (szlocal-sum(newindividual))/length(fracindices);
newindividual( fracindices(rand(1,length(fracindices)) <= p) ) = 1;

if sum(newindividual)==0, error('Empty mutation.'); end

end