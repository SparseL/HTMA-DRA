clc,clear
%No. of corresponding state vector squence
SV = 4;
% density
N = 0.2;
%Number of nodes
T = 40;
% number of data length
K = 6;
% generate density N, nodes = T, datalength = SV*(K).
% A = zeros(K*SV,T);W = zeros(T,T);
rand('twister',sum(100*clock));
for realization = 1:1
% [A,y,W] = data_generateFCM(SV,N,T,K);
% 
% savefile=sprintf('data-Node_%d_density_%d_SV_%d',T,N*100,SV);
% save(savefile,'W','y','A');
savefile=sprintf('data-Node_%d_density_%d_SV_%d',T,N*100,SV);
load(savefile);
%-----------------------MASTNet starting-------------------------
norm_index=0;
pop=100;MAX_interation=1200;
pool_size = 40;tour_size = 2;
tryCount=5;
fres = zeros(1,T);
for tryIndex = 1:tryCount
P = zeros(pop,T,T);
fit = zeros(pop,T);
Pini = zeros(pop,T,T);fitini = zeros(pop,T);
fit2ini = zeros(pop,T);
% the parameter for local search size
lambda = 1e-4;
xp_MA = zeros(T,T);
budget = (pop+2*pool_size*1500)*T;spendbudget = 0;
tic
for index = 1:T
% for each cycle, spend = pop+2*pool*gen
[a,b,c]=MA_Initialization(A,...
 y(:,index),norm_index,pop,lambda);
[xp_MA1(:,index)] = MA_Local_Search(A,y(:,index),...
norm_index,MAX_interation,pop,pool_size,tour_size,lambda,a,b,c,[]);
end

t=toc;
fprintf('CCMA Used time: %f\n',t);
%-------------------------MASTNet end------------------------------
savefile=sprintf('CommonCC-Node_%d_density_%d_SV_%d_tryIndex_%d',T,N*100,SV,tryIndex);
save(savefile,'W','xp_MA1');
end
end