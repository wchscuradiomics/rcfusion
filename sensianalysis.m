%% run search
% clear; clc;
% load our.mat;
% 
% tic;
% iz = 0.45;
% param.our.hasintercept = false;
% param.our.szlocal = round( iz * width(CAND) );
% param.our.radconrithreshold = 05;
% param.our.fitcoefs = [0.90 0.08 0.02];
% [x2, util2] = ournomogram(CAND,CLN,li,param);
% % run('ournomogram.m'); % param.our.hasintercept = true;
% elapsed_time = toc;
% 
% matname = ['our-' num2str(iz*100,'%02d') '_' num2str(param.our.radconrithreshold,'%02d') '_' ...
%   num2str(param.our.fitcoefs(1)*100,'%02d') '_' ...
%   num2str(param.our.fitcoefs(2)*100,'%02d') '_' num2str(param.our.fitcoefs(3)*100,'%02d') '.mat'];

%% use searched subset
% % step 1/5: feature selection based on newly searched subset for radiomic features
% param.our.hasintercept = true;
% LCL = CAND(:,x2==1);
% [lclsubset,B,fitinfo,bstlindex] = fs4lcl(LCL(li.trnindices,:),li.trnlabels,li.cvcvp,param.our.alpha);
% coefs = B(:,bstlindex); intercept = fitinfo.Intercept(bstlindex);
% radscores = LCL*coefs + intercept;
% [radscores,param.radiomics.radscoremus,param.radiomics.radscoresigmas] = stdfeatures(radscores,li.trnindices);
% 
% % setp 2/5: variable selection for rad-score and clinical factors, train a new nomogram
% NEW = [radscores CLN];
% [newsubset,B2,fitinfo2] = fs4new(NEW(li.trnindices,:),li.trnlabels,li.cvcvp,param.our.alpha,param.our.clnoption);
% TRN = NEW(li.trnindices,newsubset);
% [nom,cvscores] = lrm(TRN,li.trnlabels,param.our.hasintercept,li.cvcvp,true);
% [mcvauc,aucs] = meancvaucs(li.cvcvp,li.trnlabels,cvscores);
% nomcoefs = nom.Coefficients.Estimate; if param.our.hasintercept, nomcoefs=nomcoefs(2:end); end
% 
% % step 3/5: save param
% % param.radiomics.candsubset = candsubset;
% % param.radiomics.imgparams = imgparams;
% % param.radiomics.validindices5radiomics = validindices5radiomics;
% % param.cliniomics.invclnvarnames = invclnvarnames;
% param.fs.lcl.B = B; param.fs.lcl.fitinfo = fitinfo;
% param.fs.lcl.bstlindex = bstlindex; param.fs.lcl.coefs = coefs; 
% param.fs.lcl.intercept = intercept; param.fs.lcl.lclsubset = lclsubset;
% param.fs.new.newsubset = newsubset;
% param.fs.new.B = B2; param.fs.new.fitinfo = fitinfo2; param.fs.new.newsubset = newsubset;
% 
% % step 4/5: test sets based on new subset (lclsubset for feature selection, newsubset for variable selection)
% TST = NEW(li.tstindices,newsubset);
% 
% LCLSCNO4H = CANDSCNO4H(:,x2==1);  
% radscoresscno4h = LCLSCNO4H*param.fs.lcl.coefs + param.fs.lcl.intercept;
% radscoresscno4h = stdfeatures(radscoresscno4h,param.radiomics.radscoremus,param.radiomics.radscoresigmas);
% NEWSCNO4H = [radscoresscno4h CLNSCNO4H];
% 
% TSTSCNO4H = NEWSCNO4H(:,param.fs.new.newsubset);
% 
% % step 5/5: predict
% tstscores = cpredict(nom,TST,param.clnames);
% trnscores = cpredict(nom,TRN,param.clnames);
% tstscoresscno4h = cpredict(nom,TSTSCNO4H,param.clnames);
% scno4hresult.nom = perfresult(li.scno4hlabels,tstscoresscno4h,param.performance.version);
% trnresult.nom = perfresult(li.trnlabels,trnscores,param.performance.version);
% valresult.nom = perfresult(li.trnlabels,cvscores,param.performance.version);
% tstresult.nom = perfresult(li.tstlabels,tstscores,param.performance.version);
% scno4hresult.nom = perfresult(li.scno4hlabels,tstscoresscno4h,param.performance.version);
% 
% % display result
% s = [num2str(trnresult.nom.auc(1)) 9 num2str(valresult.nom.auc(1)) 9 num2str(tstresult.nom.auc(1)) 9 ...
%   num2str(scno4hresult.nom.auc(1))];
% disp(s);
% 
% save(matname);

%% epsilon
% clear; clc;
% 
% epsilons = [-15 -10 -5 00 03 05 07 10 15];
% tstaucs = zeros(size(epsilons));
% scno4haucs = zeros(size(epsilons));
% % mats = cell(1, length(epsilons));
% for i=1:length(epsilons)
%   % mats{i} = ['epsilon/our_' num2str(epsilons(i),'%02d') '_90_08_02.mat'];
%   mat = ['epsilon\our_45_' num2str(epsilons(i),'%02d') '_90_08_02.mat'];
%   tstresult = getfield(load(mat,'tstresult'),'tstresult');
%   scno4hresult = getfield(load(mat,'scno4hresult'),'scno4hresult');
%   tstaucs(i) = tstresult.nom.auc(1);
%   scno4haucs(i) = scno4hresult.nom.auc(1);
% end
% 
% plot(epsilons/10,tstaucs,'LineWidth',1.5,'Color','blue')
% set(gcf, 'Position', [100, 100, 260, 150])
% ylim([0.7 0.9]);
% hold on;
% plot(epsilons/10,scno4haucs,'LineWidth',1.0,'Color','green');
% lgd = legend({'Independent test','External test'},'FontSize',10,Position=[0.36 0.36,0.3,0.15]);
% set(lgd, 'Box', 'off')
% lgd.Position(2) = lgd.Position(2)-0.03;
% lgd.Position(4) = lgd.Position(4)*1.6;
% box on;
% xlabel('relative offset (\epsilon,%)');
% ylabel('AUC','Position',[-1.65 0.8]);
% yticks(0.75:0.05:0.9)
% xticks(-1.5:0.5:1.5);

%% z
% clear; clc;
% 
% zs = [40 42 42.5 43.5 44 45 46 46.5 47.5 50];
% tstaucs = zeros(size(zs));
% scno4haucs = zeros(size(zs));
% % mats = cell(1, length(epsilons));
% for i=1:length(zs)
%   % mats{i} = ['epsilon/our_' num2str(epsilons(i),'%02d') '_90_08_02.mat'];
%   mat = ['z\our_' num2str(zs(i)) '_05_90_08_02.mat'];
%   tstresult = getfield(load(mat,'tstresult'),'tstresult');
%   scno4hresult = getfield(load(mat,'scno4hresult'),'scno4hresult');
%   tstaucs(i) = tstresult.nom.auc(1);
%   scno4haucs(i) = scno4hresult.nom.auc(1);
% end
% 
% plot((zs-45),tstaucs,'LineWidth',1.5,'Color','blue')
% set(gcf, 'Position', [100, 100, 260, 150])
% ylim([0.7 0.9]);
% hold on;
% plot((zs-45),scno4haucs,'LineWidth',1.0,'Color','green');
% lgd = legend({'Independent test','External test'},'FontSize',10,Position=[0.36 0.26,0.3,0.16]);
% set(lgd, 'Box', 'off')
% % lgd.Position(2) = lgd.Position(2)-0.00;
% % lgd.Position(4) = lgd.Position(4)*1.6;
% box on;
% xlabel('relative offset (z, %)');
% ylabel('AUC','Position',[-4 0.86]);
% yticks(0.75:0.05:0.9)
% xticks(-5:2.5:5);

%% weights
% clear; clc;
% 
% weights = [01 02 03 04 05 06 07 08 09];
% tstaucs = zeros(size(weights));
% scno4haucs = zeros(size(weights));
% % mats = cell(1, length(epsilons));
% for i=1:length(weights)
%   % mats{i} = ['epsilon/our_' num2str(epsilons(i),'%02d') '_90_08_02.mat'];
%   mat = ['weights\our_45_05_90_' num2str(weights(i),'%02d') '_' num2str(10-weights(i),'%02d') '.mat'];
%   tstresult = getfield(load(mat,'tstresult'),'tstresult');
%   scno4hresult = getfield(load(mat,'scno4hresult'),'scno4hresult');
%   tstaucs(i) = tstresult.nom.auc(1);
%   scno4haucs(i) = scno4hresult.nom.auc(1);
% end
% 
% plot(weights,tstaucs,'LineWidth',1.5,'Color','blue')
% set(gcf, 'Position', [100, 100, 260, 150])
% ylim([0.7 0.9]);
% hold on;
% plot(weights,scno4haucs,'LineWidth',1.0,'Color','green');
% lgd = legend({'Independent test','External test'},'FontSize',10,Position=[0.36 0.30,0.3,0.15]);
% set(lgd, 'Box', 'off')
% % lgd.Position(2) = lgd.Position(2)-0.03;
% % lgd.Position(4) = lgd.Position(4)*1.5;
% box on;
% xlabel('w_2 (%, w_1 = 10 - w_2)');
% ylabel('AUC','Position',[1 0.86]);
% yticks(0.75:0.05:0.9)
% xticks(0:2.5:10);

%% table
% clear; clc;
% mats = dir('sensitive analyses\**\*.mat');
% for i = 1:height(mats)
%   path = fullfile(mats(i).folder,mats(i).name);
%   params = strsplit(path,'_');
%   tstresult = getfield(load(path,'tstresult'),'tstresult');
%   scno4hresult = getfield(load(path,'scno4hresult'),'scno4hresult');
% 
%   z = str2double(params{2});
%   params{2} = num2str(z - 45);
% 
%   epsilon = str2double(params{3});
%   params{3} = num2str(epsilon/10);
% 
%   params{6} = replace(params{6},'.mat','');
% 
%   s = [params{2} 9 params{3} 9 params{4} 9 params{5} 9 params{6} 9 ...
%     num2str(tstresult.nom.auc(1),'%4.3f') 9 num2str(scno4hresult.nom.auc(1),'%4.3f')];
%   disp(s);
% end