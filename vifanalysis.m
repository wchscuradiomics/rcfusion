clear; clc;
fsnames = {'swr.mat','stdga.mat','relieff.mat','mrmr.mat','lasso.mat','inf.mat','il.mat','rfe.mat','our.mat'}; % 
format = '%4.3f';

for i=1:length(fsnames)
  TRN = getfield(load(fsnames{i},'TRN'),'TRN');
  LCL = getfield(load(fsnames{i},'LCL'),'LCL');
  viflcl = vif(LCL);
  vifrc = vif(TRN);
  disp([fsnames{i} 9 sprintf('%.2f ', viflcl) 13 ...
   num2str(length(viflcl)) 9 num2str(sum(viflcl<2.5)) 9 num2str(sum(viflcl<5)) 13 ...
   sprintf('%.2f ', vifrc) 13 ...
   num2str(length(vifrc)) 9 num2str(sum(vifrc<2.5)) 9 num2str(sum(vifrc<5))]);
end

%%

% X = TRN;
% [n, p] = size(X);
% VIF = zeros(p, 1);
% 
% for i = 1:p
%     % 将第 i 个变量作为因变量，其余变量作为自变量
%     y_i = X(:, i);
%     X_others = X(:, [1:i-1, i+1:end]);
% 
%     % 添加常数项并进行多元线性回归
%     X_with_const = [ones(n, 1), X_others];
%     b = regress(y_i, X_with_const); 
% 
%     % 计算决定系数 R^2
%     y_hat = X_with_const * b;
%     SS_res = sum((y_i - y_hat).^2);
%     SS_tot = sum((y_i - mean(y_i)).^2);
%     R2 = 1 - (SS_res / SS_tot);
% 
%     % 计算 VIF
%     VIF(i) = 1 / (1 - R2);
% end
% disp(VIF);