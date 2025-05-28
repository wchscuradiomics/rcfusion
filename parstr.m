function value=parstr(s,index,delimiter)

if nargin==1
  index = 1;
  delimiter = '|'; 
elseif nargin==2
  delimiter = '|';
end

arr = split(s,delimiter);
value = arr{index};