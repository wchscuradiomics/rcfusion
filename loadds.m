%% load clinical datasets and basic settings
clear;clc;
load('ich3mmRS.mat'); % Center 1; training and test datasets in variable li
param.imgparamnames = {'kvp'};
imgparams = cliniomics(:,param.imgparamnames);
CLN = table2array(cliniomics);
[CLN,param.cliniomics.zscoremus,param.cliniomics.zscoresigmas] = stdfeatures(CLN,li.trnindices);
[CLN,param.cliniomics.invalidFeatureIndices] = removeNanInfFeatures(CLN);
if ~isempty(param.cliniomics.invalidFeatureIndices), error('There is empty or nan values in clinical variables.'); end
param.clnames = [1 0];
param.performance.version = 'optimal:OPTROCPT|nboot:100';

%% independent test set
wchscucliniomics = cliniomics(li.tstindices,:);
CSV = readtable(param.csv);
CSV = CSV(li.tstindices,:);
wchscuhimgparams = wchscucliniomics(:,param.imgparamnames);
RADWCHSCU = ourpreprocess(CSV,wchscuimgparams,param.radiomics);
wchsculabels = wchscucliniomics.outcome3m;
wchscucliniomics(:,param.cliniomics.invclnvarnames) = [];
CLNWCHSCU = table2array(wchscucliniomics);
CLNWCHSCU = stdfeatures(CLNWCHSCU,param.cliniomics.zscoremus,param.cliniomics.zscoresigmas);
if ~isempty(param.cliniomics.invalidFeatureIndices), CLNWCHSCU(:,param.cliniomics.invalidFeatureIndices)=[]; end
for j=1:width(CLNWCHSCU) % fill nan values
  if any(isnan(CLNWCHSCU(:,j)))
    CLNWCHSCU(isnan(CLNWCHSCU(:,j)),j) = mean(CLN(li.trnindices,j));
  end
end

%% external test set
scno4himgparams = scno4hcliniomics(:,param.imgparamnames);
RADSCNO4H = ourpreprocess(param.csv,wchscuimgparams,param.radiomics);
scno4hlabels = scno4hcliniomics.outcome3m;
scno4hcliniomics(:,param.cliniomics.invclnvarnames) = [];
CLNSCNO4H = table2array(scno4hcliniomics);
CLNSCNO4H = stdfeatures(CLNSCNO4H,param.cliniomics.zscoremus,param.cliniomics.zscoresigmas);
if ~isempty(param.cliniomics.invalidFeatureIndices), CLNSCNO4H(:,param.cliniomics.invalidFeatureIndices)=[]; end
for j=1:width(CLNSCNO4H) % fill nan values
  if any(isnan(CLNSCNO4H(:,j)))
    CLNSCNO4H(isnan(CLNSCNO4H(:,j)),j) = mean(CLN(:,j));
  end
end