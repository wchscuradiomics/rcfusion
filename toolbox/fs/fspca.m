function [fcn,TRN,TST]=fspca(XTRN,ev,XTST)
%[fcn,TRN,TST]=fspca(XTRN,ev,XTST) reduces dimension with PCA.
%
% ev: a double in (0,1) representing explained variance to keep as fraction. For example, ev = 0.95/0.98/0.92.
%
% fcn: PCA transformation function. For example, fcn(D)

if nargin == 1  
  ev = 95/100;
  XTST = [];
elseif nargin == 2
  XTST = [];
end

[pcacoefficients, pcascores, ~, ~, explained, pcacenters] = pca(XTRN);
compcount2keep = find(cumsum(explained)/sum(explained) >= ev, 1);
pcacoefficients = pcacoefficients(:,1:compcount2keep);

TRN =  pcascores(:,1:compcount2keep);
fcn = @(x)  (x - pcacenters) * pcacoefficients;
if isequal(fcn(XTRN),TRN)
  error('Not correct for pcaTransformationFcn.');
end

if ~isempty(XTST)
  TST = fcn(XTST);
else
  TST = [];
end