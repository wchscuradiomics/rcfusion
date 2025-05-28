function r = perfresult(labels,SC,option)
%r=perfresult(labels,SC,clnames,option) calculates performance result for a binary classification. 1 represents the positive
% class and 0 represents the negative class. SC(:,1) is the scores of positive samples.
%
% labels: a n-by-1 vector representing labels of n samples. 1 represents the positive class and 0 represents the negative class.
%
% SC: a n-by-2 matrix of double type, representing scores of n samples with 2 classes. SC(:,1) is the scores of positive samples.
%
% option: a string representing name-value pairs, for example, 'optimal:youden|nboot:100', 'optimal:none|nboot:0',
% 'optimal:threshold0.6|nboot:100', 'nboot:0'. 
%
% r: a structure represening the values of indicators including accuracy, sensitivity, specificity, discriminative power, tprates,
% fprates, AUC, and so on.

% setp 1: check parameters
clnames = [1 0];
if ~isequal(unique(labels),[0;1]), error('labels cannot contain values other than 0 and 1'); end
r.SCORE = SC; r.labels = labels;
bootype = 'norm'; % bca: default; or norm

% step 2: parse options
options = strsplit(option,'|');
optnames = cell(1,length(options)); optvalues = cell(1,length(options));
for i=1:length(options)
  namevalues = strsplit(options{i},':');
  optnames{i} = namevalues{1}; optvalues{i} = namevalues{2};
end

% step 3: check how to find an optimal point (for accuracy, sensitivity, specificity, etc.): none, youden, OPTROCPT, or
% accthreshold. Currently, optpoint is 'none', 'youden', 'OPTROCPT', or 'thresholdXXX'. For 'OPTROCPT', the default calculation
% equals to youden and set parameter Cost can change its calculation. If accuracy, sensitivity, specificity, etc., are not needed,
% please set '...|optimal:none|...' or do not contain 'optimal' in the parameter 'version'.  Otherwise, must find the optimal point
% and there are three ways: youden, OPTROCPT, and threshold. Refer to: https://www.mathworks.com/help/stats/perfcurve.html.
optimalindex = find(contains(optnames,'optimal','IgnoreCase',true));
if isempty(optimalindex), optpoint='none'; else, optpoint=optvalues{optimalindex}; end

% check how to calculate confidence interval
bootindex = find(contains(optnames,'nboot','IgnoreCase',true));
if isempty(bootindex) || strcmpi(optvalues{bootindex},'none')
  nbootvalue = 0; % nbootvalue=0 means no confidence interval is needed.
else
  nbootvalue = str2double(optvalues{bootindex});
end

% calculate performance indicators
if ~contains(optpoint,'threshold')
  [fprates,tprates,thresholds,auc,OPTROCPT] = perfcurve(labels,SC(:,1),clnames(1),'NBoot',nbootvalue,...
    'BootType',bootype,'Options',statset('UseSubstreams',true,'Streams',RandStream.create('mlfg6331_64')));
else % for threshold provided
  [~,~,thresholds] = perfcurve(labels,SC(:,1),clnames(1),'NBoot',nbootvalue,...
    'BootType',bootype,'Options',statset('UseSubstreams',true,'Streams',RandStream.create('mlfg6331_64')));
  threshold = str2double(replace(optpoint,'threshold',''));
  thresholds = [thresholds; threshold];
  [fprates,tprates,thresholds,auc,OPTROCPT] = perfcurve(labels,SC(:,1),clnames(1),'NBoot',nbootvalue,...
    'TVals',thresholds, ...
    'BootType',bootype,'Options',statset('UseSubstreams',true,'Streams',RandStream.create('mlfg6331_64')));
end

%%%%%%%%%% calculate accuracy based on different optimal ways %%%%%%%%%%
if contains(optpoint,'none')
  r.fprates = fprates;
  r.tprates = tprates;
  r.auc = auc;
  r.threshold = [];
  return;
elseif contains(optpoint,'OPTROCPT')  % the default of OPTROCPT is youden index, set the parameter Cost to change it.
  optrocptindex = find(fprates(:,1) == OPTROCPT(1) & tprates(:,1) == OPTROCPT(2)); % The variable OPTROCPT refers to line 80.
  optrocptindex = optrocptindex(1); % convert optrocptindex to a scalar (maybe multiple optimal ROC points exist)
  threshold = thresholds(optrocptindex,1);
elseif contains(optpoint,'youden')
  y = tprates + 1 - fprates - 1;
  [~,optrocptindex] = max(y(:,1));
  optrocptindex = optrocptindex(1);
  threshold = thresholds(optrocptindex,1);
elseif contains(optpoint,'threshold')
  optrocptindex = find(thresholds==threshold);
  optrocptindex = optrocptindex(1);
else
  error(['Not supported optimal way: ' optpoint]);
end

if contains(optpoint,'OPTROCPT') || contains(optpoint,'youden') || contains(optpoint,'threshold')
  if nbootvalue > 0 % calculate accuracy
    if contains(optpoint,'OPTROCPT') || contains(optpoint,'youden')
      [accs,~] = perfcurve(labels,SC(:,1),clnames(1),'NBoot',nbootvalue,'XCrit','accu', ...
        'BootType',bootype,'Options',statset('UseSubstreams',true,'Streams',RandStream.create('mlfg6331_64')));
    elseif contains(optpoint,'threshold')
      [accs,~] = perfcurve(labels,SC(:,1),clnames(1),'NBoot',nbootvalue,'XCrit','accu', ...
        'TVals',thresholds,...
        'BootType',bootype,'Options',statset('UseSubstreams',true,'Streams',RandStream.create('mlfg6331_64')));
    end
    accuracy = accs(optrocptindex,:);
  else
    tp = sum(labels==clnames(1) & SC(:,1)>=threshold);
    fp = sum(labels==clnames(2) & SC(:,1)>=threshold);
    npositives = sum(labels==clnames(1));
    nnegatives = sum(labels==clnames(2));
    fn = npositives - tp;
    tn = nnegatives - fp;
    accuracy = (tp+tn)/(npositives+nnegatives);
  end % end for calculating accuracy
  sensitivity = tprates(optrocptindex,:); % calculate sensitivity
  if size(fprates,2) == 3 % calculate specificity
    specificity = 1-fprates(optrocptindex,[1 3 2]);
  else
    specificity = 1-fprates(optrocptindex,:);
  end % end for calculating specificity
  dc = (sqrt(3)/pi)*(log(sensitivity./(1-sensitivity))+log(specificity./(1-specificity))); % calculate DC

  r.auc = auc;
  r.threshold = threshold;
  r.tprates = tprates;
  r.fprates = fprates;
  r.accuracy = accuracy;
  r.sensitivity = sensitivity;
  r.specificity = specificity;
  r.dc = dc;
end