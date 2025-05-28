function [newsubset,B2,fitinfo2] = fs4new(NEW,labels,cvcvp,alpha,clnoption)

if nargin==4, clnoption='lasso|mse'; end

if contains(clnoption,'lasso')
  [B2,fitinfo2] = classo(NEW,labels,'CV',cvcvp,'Alpha',alpha);
  if strcmpi(parstr(clnoption,2),'mse'),bstlindex2=fitinfo2.IndexMinMSE; else, bstlindex2=fitinfo2.Index1SE; end
  newsubset = find(B2(:,bstlindex2)~=0)';
  while length(newsubset) < 4 % Be sure clinical variables can be included. Worse performance for 1, 2, or 3.
   bstlindex2 = bstlindex2 - 1;
    newsubset = find(B2(:,bstlindex2)~=0)';
  end  
  fitinfo2.bestlambdaindex = bstlindex2;
elseif contains(clnoption,'univaranalysis')
  [~,~,newsubset] = univaranalysis(NEW,labels,'regression',str2double(parstr(clnoption,2,'|')));
  B2 = []; fitinfo2 = [];
else  
  error('Invalid clinical option in matrix NEW.');
end
