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
rand('twister',sum(100*clock));
for realization = 1:1
% the part for generate FCMs and synthetic data
% A = zeros(K*SV,T);W = zeros(T,T);
% [A,y,W] = data_generateFCM(SV,N,T,K);
% % store the useful information
% savefile=sprintf('data-Node_%d_density_%d_SV_%d',T,N*100,SV);
% save(savefile,'W','y','A');
savefile=sprintf('data-Node_%d_density_%d_SV_%d',T,N*100,SV);
load(savefile);
% %-----------------------Predefined Parameter starting-------------------------
% the parameter for the optimizer
pop=100;MAX_interation=100;
pool_size = 50;tour_size = 2;
tryCount=5;norm_index=0;lambda = 1e-4;
budget = 1200*T;
for tryIndex = 1:1:tryCount
% the intialized variable for CC
fres = zeros(1,T);P = zeros(pop,T,T);
fit = zeros(pop,T);fit2 = zeros(pop,T);
Pini = zeros(pop,T,T);fitini = zeros(pop,T);
fit2ini = zeros(pop,T);countobj = zeros(1,T);
% additional information
spendindex = zeros(1,T);countEvaid = ones(T,1);
xp_MA = zeros(T,T);
pe = 0.0; % probability for explore
% budget = (pop+pool_size*1000)*T;
% %-----------------------the proposed CC starting-------------------------
tic
i = 1;spendbudget = 0;flag = 1;
while (spendbudget < budget)
flagr = (rand<pe);
if spendbudget == 0 || flagr || flag
flag = 0;
flags = (spendbudget == 0);
for index = 1:T
if flags
    [Pini(:,:,index),fitini(:,index),fit2ini(:,index)]=MA_Initialization(A,...
        y(:,index),norm_index,pop,lambda);
    a = Pini(:,:,index);b = fitini(:,index);
    c = fit2ini(:,index);d = [];
else
    a = P(:,:,index);b = fit(:,index);
    c = fit2(:,index);d = [];
end
% for each cycle, spend = pop+2*pool*gen P(:,:,index)
[xp_MA(:,index),P(:,:,index),plot_P_best,fit(:,index),spend,fit2(:,index)] = ...
MA_Local_Search(A,y(:,index),norm_index,600,pop,pool_size,tour_size,...
lambda,a,b,c,d);
fres_temp = (plot_P_best(1)-plot_P_best(end));
if fres_temp ~= 0
    fres(index) = fres_temp;
end
% spendbudget = spendbudget+spend;
spendbudget = spendbudget+600;
spendindex(1,index) = spendindex(1,index) + spend;
fprintf('index: %d\n',index);
end
end

if length(unique(fres)) == 1
    flag = 1;
else
    [sortedF,I] = sort(fres,'descend');
    nextobj = I(1);
end
if countobj(nextobj) >= 4
   fres(nextobj) = 0;
   nextobj = I(2);
end
i = i+1;
fprintf('selected nextobj: %d\n',nextobj);
[xp_MA(:,nextobj),P(:,:,nextobj),plot_P_best,fit(:,nextobj),spend,fit2(:,nextobj)] = ...
MA_Local_Search(A,y(:,nextobj),norm_index,MAX_interation,...
pop,pool_size,tour_size,lambda,P(:,:,nextobj),fit(:,nextobj),fit2(:,nextobj),[sortedF(1:2)]);
fres_temp = (plot_P_best(1)-plot_P_best(end));
if fres_temp ~= 0
    fres(nextobj) = fres_temp;
    countobj(nextobj) = 0;
else
    countobj(nextobj) = countobj(nextobj)+1;
end

% spendbudget = spendbudget+spend;
spendbudget = spendbudget+100;
spendindex(1,nextobj) = spendindex(1,nextobj) + spend;
countEvaid(nextobj,1) = countEvaid(nextobj,1) + 1;
end

t=toc;
savefile=sprintf('CCMA1-Node_%d_density_%d_SV_%d_tryIndex_%d',T,N*100,SV,tryIndex);
save(savefile,'W','xp_MA');
fprintf('CCMA Used time: %f\n',t);
end
%-------------------------the proposed CC end------------------------
end