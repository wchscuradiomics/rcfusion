%% load clinical datasets and basic settings
clear;clc;
run('loadds.m');

%% run abl experiment: used for GLSP (ablation experiment)
% experiment with various parameter combinations and select the optimal combination
difsigtstthresholds = [0.01 0.05 0.1]; 
corrcoefthresholds = [0.6 0.8]; 
alphas = [1]; 
times = [0.45];
csvs = dir(fullfile('xxx\WCHSCU\*','*.csv'));
csvs = fullfile({csvs(:).folder}',{csvs(:).name}');

for icsv=1:length(csvs)
  for ialpha = 1:length(alphas)
    for idifsigtstthreshold = 1:length(difsigtstthresholds)
      for icorrcoefthreshold = 1:length(corrcoefthresholds)
        for itime = 1:length(times)
          % preprocess radiomics features (imaging parameters are used as radiomics features)
          [RAD,setting] = ourpreprocess(csvs{icsv},imgparams,li);
          param.radiomics = setting.radiomics;

          % preselection to obtain candidate features
          param.radiomics.preoptions = {{'ttest','corrcoef-remove'},...
            [difsigtstthresholds(idifsigtstthreshold),corrcoefthresholds(icorrcoefthreshold)]};
          candmainsubset = fs3multistages(RAD(li.trnindices,1:(width(RAD)-length(param.radiomics.imgparamindices))),li.trnlabels,...
            param.radiomics.preoptions{1},param.radiomics.preoptions{2});
          CAND = [RAD(:,candmainsubset) RAD(:,param.radiomics.imgparamindices)]; % candidate feature matrix
          candsubset = [candmainsubset param.radiomics.imgparamindices];

          % set nomogram parameters and run ablnomogram function or script
          param.abl.alpha = alphas(ialpha);
          param.abl.preoptions = param.radiomics.preoptions;
          param.abl.hasintercept = false; % false for finding and true for the final nomogram
          param.abl.szlocal = round( times(itime) * width(CAND) ); % times(itime) * width(CAND)   
          param.abl.clnoption = 'lasso|mse'; % 'lasso|mse'
          param.abl.idfmax = 0.75 ; % 0.9 0.8 0.7 0.6
          % the same setting from our          
          param.our.alpha = param.abl.alpha;
          param.our.preoptions = param.abl.preoptions;
          param.our.hasintercept = param.abl.hasintercept;
          param.our.szlocal = param.abl.szlocal;
          param.our.clnoption = param.abl.clnoption;          
          param.seed = 20240620; % 20240620
          if param.abl.szlocal/length(candsubset) > 0.7, continue; end
          [~, util] = ablnomogram(CAND,CLN,li,param);
          % run('ablnomogram.m');
          disp([csvs{icsv} 9 num2str(alphas(ialpha)) 9 num2str(times(itime)) 9 ...
            num2str(difsigtstthresholds(idifsigtstthreshold)) 9 num2str(corrcoefthresholds(icorrcoefthreshold))]);
          disp(['------------------------------------------------------------' 13]);
        end
      end
    end
  end
end

%% final FS and training, save abl.mat
clc;
%%% run section 2 with the optimal parameter combination (change function ablnomogram to a script file)
param.radiomics.candsubset = candsubset;
param.radiomics.imgparams = imgparams;
param.cliniomics.invclnvarnames = invclnvarnames;
param.fs.lcl.B = B; param.fs.lcl.fitinfo = fitinfo;
param.fs.lcl.bstlindex = bstlindex; param.fs.lcl.coefs = coefs; 
param.fs.lcl.intercept = intercept; param.fs.lcl.lclsubset = lclsubset;
param.fs.new.newsubset = newsubset;
param.fs.new.B = B2; param.fs.new.fitinfo = fitinfo2; param.fs.new.newsubset = newsubset;
param.csv = csvs{icsv};
clearvars -except param CAND cliniomics CLN LCL li NEW nom RAD radscores TRN trnresult TST tstresult ...
  util valresult x;
% save abl.mat

%% independent test
clear; clc; load('abl.mat'); warning('off');
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
% save abl.mat

%% external test
clear; clc; load('abl.mat'); warning('off');
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
% save abl.mat