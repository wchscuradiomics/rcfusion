function stat = boothl(x,y,g) % Hosmer–Lemeshow test

[predicted,observed,pvalue,chi2] = hl(x,y,g); 
stat.predicted = predicted;
stat.observed =observed;
stat.pvalue = pvalue;
stat.chi2 = chi2;

end