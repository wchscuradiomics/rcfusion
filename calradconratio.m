function radconratio = calradconratio(radconratio,radconrithreshold)

if radconratio > radconrithreshold
  radconratio = radconrithreshold-radconratio+radconrithreshold;
end

end