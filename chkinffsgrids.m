function T=chkinffsgrids(T)
%T=chkinffsgrids(T) check the argument/parameter grid of InfFS.

rmvindices = zeros(height(T),1); k = 0;
for i=1:height(T)
  if strcmpi(T.ways(i), 'U')
    T.alpha2s(i) = 0; % alpha2 and alpha3 are not needed.
    T.alpha3s(i) = 0;
  elseif strcmpi(T.ways(i), 'S')
    if T.alpha1s(i) + T.alpha2s(i) > 1, k=k+1; rmvindices(k)=i; continue; end
    T.alpha3s(i) = 1 - T.alpha1s(i) - T.alpha2s(i);
  else
    error(['Not supported way: ' T.ways(i)])
  end
end

rmvindices = rmvindices(1:k);
T(rmvindices,:) = [];
T = unique(T,'stable');