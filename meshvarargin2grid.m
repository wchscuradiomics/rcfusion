function M=meshvarargin2grid(varargin)
%M=meshvarargin2grid(varargin) join multiple variables into a m*nargin cell,
%where m is II(length(varargin{j})) and j is from 1 to nargin.
%
% varargin: an element of varargin is a column vector.

vectors = varargin;
n = nargin;
varnames = cell(1,n);
for i=1:n, varnames{i}=inputname(i); if isa(vectors{i},'double'), vectors{i} = vectors{i}(:); end, end

lens = zeros(1,n);
for j=1:n
  if ~isa(vectors{j},'cell'), vectors{j} = num2cell(vectors{j},2);end
  vectors{j} = vectors{j}(:);
  lens(j) = length(vectors{j});
end

ind = fullfact(lens);

M = cell(size(ind,1),n);
% M(:,1) = vectors{1}(ind(:,1));
for j=1:n
  M(:,j) = vectors{j}(ind(:,j));
end

M = cell2table(M,'VariableNames',varnames);

end