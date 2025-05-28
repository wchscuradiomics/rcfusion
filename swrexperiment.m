%% load clinical datasets and basic settings
clear;clc;
run('loadds.m');

%% run SWR experiment
% experiment with various parameter combinations and select the optimal combination
difsigtstthresholds = [0.01 0.05 0.1]; 
corrcoefthresholds = [0.6 0.8]; 
alphas = [1]; 
csvs = dir(fullfile('xxx\WCHSCU\*','*.csv'));
csvs = fullfile({csvs(:).folder}',{csvs(:).name}');

nomfsmethods = {'lasso','univaranalysis|0.1','univaranalysis|0.01','univaranalysis|0.05'};
criterions = {'Deviance','bic','rsquared','aic'}; % 'adjrsquared', 'sse'
prelexptparams = ctable2struct(meshvarargin2grid(csvs,difsigtstthresholds,corrcoefthresholds,alphas,...
  criterions,nomfsmethods)); % criterias

xlsxname = 'swr_prelexptresults.xlsx';
trnresults = cell(height(prelexptparams),1);
valresults = cell(size(trnresults));
tstresults = cell(size(trnresults));
lenradscores = zeros(size(trnresults)); lennoms = zeros(size(trnresults));
parfor i=1:height(prelexptparams) 
  warning('off');
  param4par = param;
  prelexptparam = prelexptparams(i);

  % preprocess radiomics features (imaging parameters are used as radiomics features)
  [RAD,setting] = ourpreprocess(prelexptparam.csv,imgparams,li);
  param4par.radiomics = setting.radiomics;

  % preselection to obtain candidate features
  param4par.radiomics.preoptions = {{'ttest','corrcoef-remove'},...
    [prelexptparam.difsigtstthreshold,prelexptparam.corrcoefthreshold]};
  candmainsubset = fs3multistages(RAD(li.trnindices,1:(width(RAD)-length(param4par.radiomics.imgparamindices))),li.trnlabels,...
    param4par.radiomics.preoptions{1},param4par.radiomics.preoptions{2});
  CAND = [RAD(:,candmainsubset) RAD(:,param4par.radiomics.imgparamindices)]; % candidate feature matrix
  candsubset = [candmainsubset param4par.radiomics.imgparamindices];

  % set nomogram parameters and run swrnomogram function or script
  param4par.swr.alpha = prelexptparam.alpha;
  param4par.swr.hasintercept = true; % false for finding and true for the final nomogram
  param4par.swr.clnoption = 'lasso|mse'; % 'lasso|mse' 
  param4par.seed = 20240620; % 20240620
  
  [trnresult,valresult,lclsubset,newsubset,param4par] = swrnomogram(CAND,CLN,li,param4par,prelexptparam);
   
  trnresults{i} = trnresult; valresults{i} = valresult;
  lenradscores(i) = length(lclsubset); lennoms(i) = length(newsubset);
end

[trnresults,valresults,lenradscores,lennoms,prelexptparams] = rmmtycells(trnresults,...
  trnresults,valresults,lenradscores,lennoms,prelexptparams);
prelexptresults = wrprelresults(xlsxname,prelexptparams,trnresults,valresults,table(lenradscores,lennoms));


%% final FS and training, save swr.mat
clc;
%%% run section 2 with the optimal parameter combination (change function swrnomogram to a script file) 
param.radiomics.candsubset = candsubset;
param.radiomics.imgparams = imgparams;
param.cliniomics.invclnvarnames = invclnvarnames;
param.fs.lcl.coefs = coefs; param.fs.lcl.intercept = intercept;
param.fs.lcl.lclsubset = lclsubset;
param.fs.new.newsubset = newsubset;
if strcmpi(prelexptparam.nomfsmethod,'lasso'), param.fs.new.B = B2; param.fs.new.fitinfo = fitinfo2; end
param.fs.new.newsubset = newsubset;
param.csv = prelexptparam.csv;
clearvars -except param CAND cliniomics CLN LCL li NEW nom RAD radscores TRN trnresult ...
  util valresult x;
% save swr.mat

%% independent test
clear; clc; load('swr.mat'); warning('off');
clear;clc;
run('loadds.m');

CANDWCHSCU = RADWCHSCU(:,param.radiomics.candsubset);
LCLWCHSCU = CANDWCHSCU(:,x==1);  
radscoreswchscu = LCLWCHSCU*param.fs.lcl.coefs + param.fs.lcl.intercept;
radscoreswchscu = stdfeatures(radscoreswchscu,param.radiomics.radscoremus,param.radiomics.radscoresigmas);
NEWWCHSCU = [radscoresscno4h CLNWCHSCU];
TSTWCHSCU = NEWWCHSCU(:,param.fs.new.newsubset);

tstscoreswchscu = cpredict(nom,TSTWCHSCU,param.clnames);
wchscuresult.nom = perfresult(wchsculabels,tstscoreswchscu,param.performance.version);
disp(wchscuresult.nom.auc(1));
% save swr.mat

%% external test
clear; clc; load('swr.mat'); warning('off');
clear;clc;
run('loadds.m');

CANDSCNO4H = RADSCNO4H(:,param.radiomics.candsubset);
LCLSCNO4H = CANDSCNO4H(:,x==1);  
radscoresscno4h = LCLSCNO4H*param.fs.lcl.coefs + param.fs.lcl.intercept;
radscoresscno4h = stdfeatures(radscoresscno4h,param.radiomics.radscoremus,param.radiomics.radscoresigmas);
NEWSCNO4H = [radscoresscno4h CLNSCNO4H];
TSTSCNO4H = NEWSCNO4H(:,param.fs.new.newsubset);

tstscoresscno4h = cpredict(nom,TSTSCNO4H,param.clnames);
scno4hresult.nom = perfresult(scno4hlabels,tstscoresscno4h,param.performance.version);
disp(scno4hresult.nom.auc(1));
% save swr.mat