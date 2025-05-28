function varargout=rmmtycells(bcells,varargin)
%varargout=rmmtycells(bcells,varargin) removes empty elements.
%
% bcell: a cell with n-by-1 elements and it is the baseline.

vldindices = ~cellfun(@isempty,bcells);
varargout = cell(1,length(varargin));
for i=1:length(varargin)
  cells = varargin{i};
  varargout{i} = cells(vldindices,:);
end