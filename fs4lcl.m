function [lclsubset,B,fitinfo,bstlindex] = fs4lcl(LCL,labels,cvcvp,alpha,dfmax,criteria)

if nargin==4, dfmax=Inf; criteria='auto'; elseif nargin==5, criteria='auto'; end

% numtrnsamples = height(LCL);
[B,fitinfo] = lasso(LCL,labels,'CV',cvcvp,'Alpha',alpha,'DFmax',dfmax); %#ok<*PFBNS>
if strcmpi(criteria,'auto')
  if alpha<1, bstlindex=fitinfo.Index1SE; else, bstlindex=fitinfo.IndexMinMSE; end
else
  bstlindex = selbestlambda(fitinfo,criteria);
end
lclsubset = find(B(:,bstlindex)~=0)';

