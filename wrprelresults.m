function prelexptresult = wrprelresults(xlsxname,prelexptparams,trnresults,valresults,tstresults,T)

trnresults = [trnresults{:}]';
valresults = [valresults{:}]';
tstresults = [tstresults{:}]';
if nargin==5, T=[]; end

mls = fieldnames(trnresults(1));
indicators = cell(1,length(mls));
varnames = cell(1,length(mls));
for j=1:length(mls)
  jtrnresults = [trnresults(:).(mls{j})]';
  jvalresults = [valresults(:).(mls{j})]';
  jtstresults = [tstresults(:).(mls{j})]';
  indicators{j} = [[jtrnresults(:).auc]' [jtrnresults(:).accuracy]' ...
    [jvalresults(:).mcvauc]' [jvalresults(:).stdmcvauc]' [jvalresults(:).auc]' [jvalresults(:).accuracy]' ...
    [jtstresults(:).auc]' [jtstresults(:).accuracy]'];
  varnames{j} = {[mls{j} 'trnauc'],[mls{j} 'trnacc'],[mls{j} 'mcvauc'],[mls{j} 'stdmcvauc'], ...
    [mls{j} 'cvauc'] ,[mls{j} 'cvacc'], [mls{j} 'tstauc'],[mls{j} 'tstacc']};  
end
indicators = cell2mat(indicators);
varnames = horzcat(varnames{:}); 

prelexptresult = [struct2table(prelexptparams) array2table(indicators,'VariableNames',varnames)];
if ~isempty(T), prelexptresult = [prelexptresult T]; end

writetable(prelexptresult,xlsxname,'FileType','spreadsheet','Sheet','main','WriteMode','overwritesheet');