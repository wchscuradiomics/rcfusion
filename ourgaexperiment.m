%% load clinical datasets (Part 1 of loadds.m)
clear;clc;
load('ich3mmRS.mat'); % training and test datasets in variable li
cliniomics = wchscucliniomics; clear wchscucliniomics scno4hcliniomics;
validindices5radiomics = cliniomics.vldindices5pyradiomics;
param.imgparamnames = {'kvp'};
imgparams = cliniomics(:,param.imgparamnames);
cliniomics(:,invclnvarnames) = [];
CLN = table2array(cliniomics); % CLN is a matrix now
[CLN,param.cliniomics.zscoremus,param.cliniomics.zscoresigmas] = stdfeatures(CLN,li.trnindices);
[CLN,param.cliniomics.invalidFeatureIndices] = removeNanInfFeatures(CLN);
if ~isempty(param.cliniomics.invalidFeatureIndices), error('There is empty or nan values in clinical variables.'); end
param.clnames = [1 0];
param.performance.version = 'optimal:youden|nboot:100';

%% run fusion experiment
clc;
% experiment with various parameter combinations and select the optimal combination
difsigtstthresholds =  [0.05]; % difsigtstthresholds = [0.01 0.05 0.1]; 
corrcoefthresholds = [0.8]; % corrcoefthresholds = [0.6 0.8]; 
alphas = [1]; % alphas = [0.5 1];
times = [0.45];
% csvs = dir(fullfile('\\xxx\xx\**','*.csv'));
% csvs = fullfile({csvs(:).folder}',{csvs(:).name}');
% csvs = {'\\xxx\xx\ICH\3mmRS\Count32Width4WaveletWidth4Level2_rbio31_wchscu.csv'};
csvs = {'\\xxx\xx\Count16Width8WaveletWidth8Level2\coif3.csv'};

for icsv=1:length(csvs)
  for ialpha = 1:length(alphas)
    for idifsigtstthreshold = 1:length(difsigtstthresholds)
      for icorrcoefthreshold = 1:length(corrcoefthresholds)
        for itime = 1:length(times)
          % preprocess radiomics features (imaging parameters are used as radiomics features)
          [RAD,setting] = ourpreprocess(csvs{icsv},validindices5radiomics,imgparams,li);
          param.radiomics = setting.radiomics;

          % preselection to obtain candidate features
          param.radiomics.preoptions = {
            {'ttest','corrcoef-remove'},...
            [difsigtstthresholds(idifsigtstthreshold),corrcoefthresholds(icorrcoefthreshold)] ...
          };
          candmainsubset = fs3multistages(...
            RAD(li.trnindices,1:(width(RAD)-length(param.radiomics.imgparamindices))), li.trnlabels, ...
            param.radiomics.preoptions{1},param.radiomics.preoptions{2} ...
          );
          CAND = [RAD(:,candmainsubset) RAD(:,param.radiomics.imgparamindices)]; % candidate feature matrix
          candsubset = [candmainsubset param.radiomics.imgparamindices];

          % set our parameters and run our framework (nomogram script)
          param.our.alpha = alphas(ialpha);
          param.our.preoptions = param.radiomics.preoptions;
          param.our.hasintercept = false; % false for finding and true for the final nomogram
          param.our.szlocal = round( times(itime) * width(CAND) ); %  a preset local range in Eq. (4)
          param.our.clnoption = 'lasso|mse'; % 'lasso|mse'
          param.our.radconrithreshold = 05; % relative offset in Eq. (8)
          param.our.fitcoefs = [0.90 0.08 0.02]; % w_1, w_2, and w_3 in Eqs. (5)-(9)
          param.seed = 20240620;
          if param.our.szlocal/length(candsubset) > 0.7, continue; end
          [~, util] = ournomogram(CAND,CLN,li,param);         
          % run('ournomogram.m'); % param.our.hasintercept = true;
          elapsed_time = toc;
          disp([ ...
            csvs{icsv} 9 num2str(alphas(ialpha)) 9 num2str(times(itime)) 9 ...
            num2str(difsigtstthresholds(idifsigtstthreshold)) 9 num2str(corrcoefthresholds(icorrcoefthreshold)) ...
          ]);
          disp([repmat('-', 1, 72) 13]);
        end
      end
    end
  end
end

%% save our.mat
% clc;
% param.radiomics.candsubset = candsubset;
% param.radiomics.imgparams = imgparams;
% param.radiomics.validindices5radiomics = validindices5radiomics;
% param.cliniomics.invclnvarnames = invclnvarnames;
% param.fs.lcl.B = B; param.fs.lcl.fitinfo = fitinfo;
% param.fs.lcl.bstlindex = bstlindex; param.fs.lcl.coefs = coefs; 
% param.fs.lcl.intercept = intercept; param.fs.lcl.lclsubset = lclsubset;
% param.fs.new.newsubset = newsubset;
% param.fs.new.B = B2; param.fs.new.fitinfo = fitinfo2; param.fs.new.newsubset = newsubset;
% param.csv = csvs{icsv};
% clearvars -except param CAND cliniomics CLN LCL li NEW nom RAD radscores TRN trnresult TST tstresult util valresult x;
% % save our.mat

%% external validation (Part 3 of loadds.m)
% clear; clc; 
% load('our.mat'); 
% scno4hcliniomics = getfield(load('ich3mmRS.mat','scno4hcliniomics'),'scno4hcliniomics');
% validindices5radiomics = scno4hcliniomics.vldindices5pyradiomics;
% imgparams = scno4hcliniomics(:,param.imgparamnames);
% RADSCNO4H = ourpreprocess(replace(param.csv,'WCHSCU','SCNO.4H'),validindices5radiomics,imgparams,param.radiomics);
% 
% scno4hlabels = scno4hcliniomics.outcome3m;
% scno4hcliniomics(:,param.cliniomics.invclnvarnames) = [];
% CLNSCNO4H = table2array(scno4hcliniomics);
% CLNSCNO4H = stdfeatures(CLNSCNO4H,param.cliniomics.zscoremus,param.cliniomics.zscoresigmas);
% if ~isempty(param.cliniomics.invalidFeatureIndices), CLNSCNO4H(:,param.cliniomics.invalidFeatureIndices)=[]; end
% for j=1:width(CLNSCNO4H) % fill nan values
%   if any(isnan(CLNSCNO4H(:,j)))
%     CLNSCNO4H(isnan(CLNSCNO4H(:,j)),j) = mean(CLN(:,j));
%   end
% end
% 
% CANDSCNO4H = RADSCNO4H(:,param.radiomics.candsubset);
% LCLSCNO4H = CANDSCNO4H(:,x==1);  
% radscoresscno4h = LCLSCNO4H*param.fs.lcl.coefs + param.fs.lcl.intercept;
% radscoresscno4h = stdfeatures(radscoresscno4h,param.radiomics.radscoremus,param.radiomics.radscoresigmas);
% NEWSCNO4H = [radscoresscno4h CLNSCNO4H];
% TSTSCNO4H = NEWSCNO4H(:,param.fs.new.newsubset);
% 
% tstscoresscno4h = cpredict(nom,TSTSCNO4H,param.clnames);
% scno4hresult.nom = perfresult(scno4hlabels,tstscoresscno4h,param.performance.version);
% disp(scno4hresult.nom.auc(1));
% li.scno4hlabels = scno4hlabels;
% clear validindices5radiomics imgparams j noneligindices scno4hlabels ans;
% % save our.mat