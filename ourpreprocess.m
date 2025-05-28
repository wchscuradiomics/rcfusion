function [RAD,param]=ourpreprocess(csv,imgparams,li)

omics = li;

% extract radiomics features from a csv file -> RAD, a numerical table
if ~isa(csv,'double')
  RAD = readtable(csv); % RAD will be a matrix (converted from the table output.csv)
else
  RAD = csv;
end
[startindex,param.radiomics.names]=fid1var5pyradiomics(RAD.Properties.VariableNames);
RAD = RAD(param.radiomics.validindices,startindex:end);

% add imaging paramters to RAD (RAD is still a numerical table)
if isfield(li,'labels')
  param.radiomics.imgparamnames = imgparams.Properties.VariableNames; % add imaging parameters
  param.radiomics.names = [param.radiomics.names param.radiomics.imgparamnames];
  param.radiomics.imgparamindices = (width(RAD)+1):(width(RAD)+length(param.radiomics.imgparamnames));
  RAD = [RAD imgparams(:,param.radiomics.imgparamnames)];
else
  RAD = [RAD imgparams(:,omics.imgparamnames)];
end


% convert RAD to matrix and standardize RAD
RAD = ctable2array(RAD);
if isfield(li,'labels')
  [RAD,param.radiomics.invalidFeatureIndices] = rmvnaninffeatures(RAD);
  [RAD,param.radiomics.zscoremus,param.radiomics.zscoresigmas] = stdfeatures(RAD,li.trnindices);
else
  if ~isempty(omics.invalidFeatureIndices), RAD(:,omics.invalidFeatureIndices)=[]; end
  [RAD,~] = stdfeatures(RAD,omics.zscoremus,omics.zscoresigmas);
end

end