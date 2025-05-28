function d = diff2(v)

ev = [v 0 0];
d = zeros(1,length(v));

for j=1:length(v) 
  d(j) = ev(j+2) - 2*ev(j+1) + ev(j);
end