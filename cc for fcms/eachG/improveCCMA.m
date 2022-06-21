clc,clear 
%No. of corresponding state vector squence
SV = 10;
%Number of nodes
T = 40;
% number of data length
K = 8;
% generate density N, nodes = T, datalength = SV*(K).
% gen = 1000:2000:10000;
gen = 1000;NN = [0.2,0.4];
PEE = 0.05:0.05:0.05;
eachG = [10:10:100];
for density = 1:2
% density
N = NN(density);
% the part for generate FCMs and synthetic data
A = zeros(K*SV,T);W = zeros(T,T);
[A,y,W] = data_generateFCM(SV,N,T,K);
% store the useful information
savefile=sprintf('data-Node_%d_density_%d_SV_%d',T,N*100,SV);
save(savefile,'W','y','A');
% load(savefile);
for re = 1:length(eachG)
% %-----------------------Predefined Parameter starting-------------------------
% the parameter for the optimizer
pop=100;MAX_interation=gen;
eachgen = eachG(re);pt = 1;repeatITER = 200;
pool_size = 50;tour_size = 2;
tryCount=10;norm_index=0;lambda = 10e-4;
% the intialized variable for CC
fres = zeros(1,T);P = zeros(pop,T,T);
fit = zeros(pop,T);fit2 = zeros(pop,T);
Pini = zeros(pop,T,T);fitini = zeros(pop,T);
fit2ini = zeros(pop,T);countobj = zeros(1,T);
% additional information
spendindex = zeros(1,T);countEvaid = ones(T,1);
xp_MA = zeros(T,T);
nextset = [];
% pe = PEE(re); % probability for explore
pe = 0.05; 
budget = MAX_interation*T;
% %-----------------------the proposed CC starting-------------------------
tic
for tryIndex = 1:tryCount
fprintf('tryIndex: %d\n',tryIndex);
spendbudget = 0;flag = 1;
while (spendbudget < budget)
flagr = (rand<pe);
if spendbudget == 0 || flag || flagr
flag = 0;
flags = (spendbudget == 0);
% for index = 1:T
index = 1;
while (spendbudget < budget) && (index <= T)
if flags
    [Pini(:,:,index),fitini(:,index),fit2ini(:,index)]=MA_Initialization(A,...
        y(:,index),norm_index,pop,lambda);
%     [a,b,c]=MA_Initialization(A,y(:,index),norm_index,pop,lambda);
    a = Pini(:,:,index);b = fitini(:,index);
    c = fit2ini(:,index);
    d = [];
else
    a = P(:,:,index);b = fit(:,index);
    c = fit2(:,index);d = [];
end
% for each cycle, spend = pop+2*pool*gen P(:,:,index)
if countobj(index)*eachgen>=repeatITER
    index = index+1;
    continue;
else
[xp_MA(:,index),P(:,:,index),plot_P_best,fit(:,index),spend,fit2(:,index)] = ...
MA_Local_Search(A,y(:,index),norm_index,eachgen,pop,pool_size,tour_size,...
lambda,a,b,c,d);
fres_temp = (plot_P_best(1)-plot_P_best(end));
if fres_temp ~= 0
    fres(index) = fres_temp;
    countobj(index) = 0;
else
    countobj(index) = countobj(index)+1;
end
plot_P_best = [];
% spendbudget = spendbudget+spend;
spendbudget = spendbudget+eachgen;
spendindex(1,index) = spendindex(1,index) + spend;
countEvaid(index,1) = countEvaid(index,1) + 1;
% fprintf('index: %d\n',index);
end
index = index+1;
end
end
flage = 1;
while (spendbudget < budget) && flage
if length(unique(fres)) == 1
    flag = 1;flage = 0;
else
    [sortedF,I] = sort(fres,'descend');
    nextobj = I(1);
end
fprintf('selected nextobj: %d\n',nextobj);
[xp_MA(:,nextobj),P(:,:,nextobj),plot_P_best,fit(:,nextobj),spend,fit2(:,nextobj)] = ...
MA_Local_Search(A,y(:,nextobj),norm_index,eachgen,...
pop,pool_size,tour_size,lambda,P(:,:,nextobj),fit(:,nextobj),fit2(:,nextobj),[sortedF(1:2)]);
fres_temp = (plot_P_best(1)-plot_P_best(end));
if fres_temp < sortedF(2)
    flage = 0;
end
if fres_temp ~= 0
    fres(nextobj) = fres_temp;
    countobj(nextobj) = 0;
else
    countobj(nextobj) = countobj(nextobj)+1;
end
if countobj(nextobj)*eachgen >= repeatITER
   fres(nextobj) = 0;
   nextobj = I(2);
   flage = 0;
end
plot_P_best = [];
% spendbudget = spendbudget+spend;
spendbudget = spendbudget+eachgen;
spendindex(1,nextobj) = spendindex(1,nextobj) + spend;
countEvaid(nextobj,1) = countEvaid(nextobj,1) + 1;
end
end
t=toc;
savefile=sprintf('CCMA-Node_%d_density_%d_SV_%d_try_%d_re_%d',T,N*100,SV,tryIndex,re);
save(savefile,'W','xp_MA','countEvaid');
% fprintf('CCMA Used time: %f\n',t);

%-------------------------the proposed CC end-----------------------
%-------------------------Common CC starting------------------------
% xp_MA1 = zeros(T,T);
% commonbest = zeros(T,gen);
% tic
% for index = 1:T
% % for each cycle, spend = pop+2*pool*gen
% [xp_MA1(:,index),~,commonbest(index,:)] = MA_Local_Search(A,y(:,index),...
% norm_index,gen,pop,pool_size,tour_size,lambda,Pini(:,:,index),...
% fitini(:,index),fit2ini(:,index),[]);
% fprintf('commonCC_index: %d\n',index);
% end
% toc
% savefile=sprintf('CommonCC-Node_%d_density_%d_SV_%d_try_%d_re_%d',T,N*100,SV,tryIndex,re);
% save(savefile,'xp_MA1','commonbest');
%-------------------------Common CC end-----------------------------
end
end
end