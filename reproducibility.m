clear; clc;


load our.mat; % or use the following two steps:


% % step 1/2: preprocess clinical information, global settings 
% load('ich3mmRS.mat'); % training and test datasets in variable li
% cliniomics = wchscucliniomics; clear wchscucliniomics scno4hcliniomics;
% validindices5radiomics = cliniomics.vldindices5pyradiomics;
% param.imgparamnames = {'kvp'};
% imgparams = cliniomics(:,param.imgparamnames);
% cliniomics(:,invclnvarnames) = [];
% CLN = table2array(cliniomics); % CLN is a matrix now
% [CLN,param.cliniomics.zscoremus,param.cliniomics.zscoresigmas] = stdfeatures(CLN,li.trnindices);
% [CLN,param.cliniomics.invalidFeatureIndices] = removeNanInfFeatures(CLN);
% if ~isempty(param.cliniomics.invalidFeatureIndices), error('There is empty or nan values in clinical variables.'); end
% param.clnames = [1 0];
% param.performance.version = 'optimal:youden|nboot:100';
% 
% % step 2/2: preprocess radiomics information
% csv = {'xx\coif3.csv'};
% [RAD,setting] = ourpreprocess(csv,validindices5radiomics,imgparams,li);
% param.radiomics = setting.radiomics;
% 
% % preselection to obtain candidate features
% param.radiomics.preoptions = {{'ttest','corrcoef-remove'},[0.05,0.8]};
% candmainsubset = fs3multistages(RAD(li.trnindices,1:(width(RAD)-length(param.radiomics.imgparamindices))), ...
%   li.trnlabels,param.radiomics.preoptions{1},param.radiomics.preoptions{2});
% CAND = [RAD(:,candmainsubset) RAD(:,param.radiomics.imgparamindices)]; % candidate feature matrix
% candsubset = [candmainsubset param.radiomics.imgparamindices];


% set our parameters and run our framework (nomogram function)
tic;
param.our.alpha = 1;
param.our.preoptions = param.radiomics.preoptions;
param.our.clnoption = 'lasso|mse';
param.seed = 20240620;
iz = 0.45; % local search range, a ratio of width(CAND)
param.our.hasintercept = false; % set to false for convenient contribution calculation; set to true when constructing nomogram.
param.our.szlocal = round( iz * width(CAND) );
param.our.radconrithreshold = 05; % \epsilon, baseline contribution of rad-score = 61.2%
param.our.fitcoefs = [0.90 0.08 0.02]; % w1, w2, and w3
[x, util] = ournomogram(CAND,CLN,li,param);
elapsed_time = toc;