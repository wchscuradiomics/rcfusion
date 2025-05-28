function v=idxvalue(obj,index)

if isa(obj,'cell')
  v = obj{index};
elseif isa(obj, 'double')
  v = obj(index);
else
  error('Just cell or double is supported.');
end