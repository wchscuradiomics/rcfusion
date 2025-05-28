function [p,h,stats] = mwu(FP,FN,side)

if nargin == 2
  side = 'both';
end

if size(FP,2) ~= size(FN,2)
  error 'column number is not equal'
end
p = zeros(1,size(FP,2));
h = zeros(1,size(FP,2));
stats = cell(1,size(FP,2));
for i=1:size(FP,2)
  [p(i),h(i),stats{i}] = ranksum(FP(:,i),FN(:,i),'tail',side);
end