clc,clear 
%No. of corresponding state vector squence
SV = 10;
% density
N = 0.2;
%Number of nodes
T = 40;
% number of data length
K = 6;
% generate density N, nodes = T, datalength = SV*(K).
A = zeros(K*SV,T);W = zeros(T,T);
for realization = 1:1
% [A,y,W] = data_generateFCM(SV,N,T,K);

% savefile=sprintf('data-Node_%d_density_%d_SV_%d',T,N*100,SV);
% save(savefile,'W','y','A');
savefile=sprintf('data-Node_%d_density_%d_SV_%d',T,N*100,SV);
load(savefile);
% %-----------------------MASTNet starting-------------------------
norm_index=1;
pop=1;MAX_interation=500;tryCount=1;
% the parameter for local search size
xp_MA = zeros(T,T);
% budget = (pop+MAX_interation)*T;spendbudget = 0;
tic
for tryIndex = 1:tryCount
xp_MA = zeros(T,T);
for index = 1:1
% for each cycle, spend = pop+2*pool*gen
[xp_MA(:,index),P,plot_P_best,fit] = ...
fastMA_Local_Search(A,y(:,index),norm_index,MAX_interation,pop);
fprintf('index: %d\n',index);
end
toc
savefile=sprintf('Common-Node_%d_density_%d_SV_%d',T,N*100,SV);
save(savefile,'W','xp_MA');
end
end