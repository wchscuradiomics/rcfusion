function POP=ourcreation(popsize,RAD,CLN,labels,cvcvp,util)
%POP=ourcreation creates an initial population.
%
% POP: a popsize-by-nvars matrix specifying the created initialized population, where popsize is options.PopulationSize.

% if exist('INIT.txt','file'), POP = readmatrix('INIT.txt'); return; end

nvars = width(RAD); POP = nan(popsize,nvars); szlocal = util.param.our.szlocal;

for i=1:popsize
  while true
    individual = double(rand(1,nvars) <= szlocal/nvars);
    if ~chkrepeatability(POP,i-1,individual,RAD,CLN,labels,cvcvp,util)
      continue;
    else
      POP(i,:) = individual;
      break;
    end
  end
end % end for all new individuals
util.INIT = POP;
end % end for function ourcreation

%% eligible=chkrepeatability(POP,k,individual)
function eligible = chkrepeatability(POP,k,individual,RAD,CLN,labels,cvcvp,util)

eligible = true; alpha = util.param.our.alpha; hasintercept = util.param.our.hasintercept; clnoption = util.param.our.clnoption;

% check 1
if sum(individual) == 0, eligible=false; return; end

% check 2
for i=1:k
  if sum(individual==1 & POP(i,:)==individual)/sum(POP(i,:)==1) > 0.7
    eligible = false;
    return;
  end
end

% check 3
rad1indices = find(individual==1); % indices of gene = 1 in radiomics features
LCL = RAD(:,rad1indices); nlclvars = width(LCL);
[lclsubset,B,fitinfo,bstlindex] = fs4lcl(LCL,labels,cvcvp,alpha);
coefs = B(:,bstlindex); intercept = fitinfo.Intercept(bstlindex);
radscores = LCL*coefs + intercept;
radscores = stdfeatures(radscores,(1:height(radscores))');
NEW = [radscores CLN];
[newsubset,B2,fitinfo2] = fs4new(NEW,labels,cvcvp,alpha,clnoption);
if length(lclsubset)<4 || newsubset(1)~=1, eligible=false; return; end
[nom,cvscores] = lrm(NEW(:,newsubset),labels,hasintercept,cvcvp,true);
[mcvauc,~] = meancvaucs(cvcvp,labels,cvscores);
nomcoefs = nom.Coefficients.Estimate; if hasintercept, nomcoefs=nomcoefs(2:end); end
if abs(nomcoefs(1))/sum(abs(nomcoefs)) < 0.2 || mcvauc < 0.8, eligible=false; end %  || vc > 0.15
end