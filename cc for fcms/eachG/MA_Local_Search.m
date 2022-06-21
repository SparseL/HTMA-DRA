%% editer: Kai Wu
%  opitimizaton: min ||x||_1 && ||AX-y||_2^2
%  x is a sparse singal. y is the observation. A is observation matrix.
%  MOEA+Local Search
%% main function 
function [P_best,P,plot_P_best,fit,spendbudget,fit2]...
    =MA_Local_Search(A,y,norm_index,MAX_interation,popsize,pool_size,tour_size,lambda,Pstore,fitstore,fit2store,C) 
%% initialization
%-------------------initialization start-----------------------------------

[~,N]=size(A);gen = 1;
spendbudget = 0;randlocal = 0.2;
P = Pstore;
fit = fitstore;
fit2 = fit2store;
min_ = find(fit == min(fit));
P_best = P(min_(1),:);
plot_P_best(gen) = fit(min_(1));
P_best_fit2 = fit2(min_(1));
flag = 1;
%-------------------initialization end-------------------------------------
%% main loop  (gen<=MAX_interation) && (flag) && (P_best_fit2 >=1e-6)
%%%-----------------main loop start----------------------------------------
while (gen<=MAX_interation) && (flag)
% the parameter needs to be optimized. 
% fprintf('CCMA generations: %d\n',gen);
%% selection   
    P_selection = tournament_selection(P, pool_size, tour_size, fit);
    %% crossover
    crossover_rate=0.5;
    sub_cross=ceil(crossover_rate*pool_size);
    if mod(sub_cross,2)==1
        sub_cross=sub_cross+1;
    end
     [~,sortcross]=sort(rand(1,pool_size));
     sub_crossP=P_selection(sortcross(1:1:sub_cross),:);
     for cr=1:1:sub_cross/2     
         %-----------------BLX-alpha---------------------------------------
         ALPHA_crossover=rand(1,N);
         differ=sub_crossP(cr,:)-sub_crossP(sub_cross-cr,:);
         di=abs(differ);min_differ=min([sub_crossP(cr,:);sub_crossP(cr,:)]);
         X_i_1=min_differ-ALPHA_crossover.*di;
         X_i_2=min_differ+ALPHA_crossover.*di;
         sub_crossP(cr,:)=X_i_1;
         sub_crossP(sub_cross-cr,:)=X_i_2;
         clear X_i_2 X_i_1;
         %-----------------BLX-alpha---------------------------------------
     end
    sub_crossP = [sub_crossP;P_selection(sortcross(sub_cross+1:1:end),:)];
    %% mutation
    mutation_rate=0.1;
    mutation_mount=ceil(pool_size*N*mutation_rate);
    mutation_sub_index=randi([1,pool_size*N],1,mutation_mount);
    muP=sub_crossP;
    ii=(1-gen/MAX_interation).^2;
    J_J=rand.^ii;
    for jj=1:1:mutation_mount
        muP(mutation_sub_index(jj))=muP(mutation_sub_index(jj))+(1-J_J);
    end
%     muP = deleteNeg(muP);
    %muP is the population
    % the part for find unique solution
    ii = 1;newi = [];
    for i2 = 1:size(muP,1)-1
        if notExist(muP(i2,:),[],P,muP,i2+1)==1        %avoid repetition
            newi(ii) = i2;
            ii = ii+1;
        end
    end
    newi(ii) = size(muP,1);
    muP_final = muP(newi,:);
    
    [fit_muP,muP_fit2]=fitness(muP_final,lambda);
    spendbudget = spendbudget + size(muP_final,1);
    %% local search   
%------^^^^^-------local search----------------------------------------
    % case one 
    newP = [];
    II = randperm(length(fit_muP));
    for i3 = 1:floor(randlocal*length(fit_muP))
         newPP = halflocalsearch(A, y, muP_final(II(i3),:));
         newP = [newP;newPP];
    end
%     newP=deleteNeg(newP);
%     newP=Local_Search(muP_final,fit_muP,y,A,lambda);
%     newP=deleteNeg(newP,index);
    newii = [];
    ii = 2;newii(1) = size(newP,1);
    for i4 = 1:1:(size(newP,1)-1)
        if notExist(newP(i4,:),muP_final,P,newP,i4+1)==1        %avoid repetition
            newii(ii) = i4;
            ii = ii+1;
        end
    end
    newP_final = newP(newii,:);  
    [newP_fit,newP_fit2] = fitness(newP_final,lambda);
    spendbudget = spendbudget + size(newP_final,1); 
%-----------------update population------------------------------------
    combinedP=[newP_final;muP_final;P];
    fit_combinedP=[newP_fit;fit_muP;fit];
    fit_combinedP2=[newP_fit2;muP_fit2;fit2];
    [~,LC]=sort(fit_combinedP,'ascend');
    P = combinedP(LC(1:popsize),:);
    fit = fit_combinedP(LC(1:popsize));
    fit2 = fit_combinedP2(LC(1:popsize));
%-----------------update best indidual-------------------------------------
    clear min_
    min_ = find(fit == min(fit));
    P_best = P(min_(1),:);
    plot_P_best(gen+1) = fit(min_(1));
    P_best_fit2 = fit2(min_(1),:);
    gen=gen+1;
%     if ~isempty(C)
%         flag = ((plot_P_best(1)-plot_P_best(end)) < C(2));
%     else
%         flag = 1;
%     end
%%%-----------------main loop end------------------------------------------
end
%% subfunction for calculate fitness
 function [fit_,fit2]=fitness(P_,lambda)
        [popsize_,~]=size(P_);
        fit1=zeros(popsize_,1);
        if norm_index==1
            fit1(:,1)=sum(abs(P_),2);
        elseif norm_index==0
            e=zeros(size(P_));
            e(P_~=0)=1;
            fit1(:,1)=sum(e,2);
        elseif norm_index==1/2
            temp=(sum(sqrt(abs(P_)),2));
            fit1(:,1)=temp.^2;
        end
        fit2=sum((abs(A*P_'-repmat(y,[1,popsize_]))).^2,1)/(2*length(y));
         fit2= fit2';
         fit_ = fit2+lambda*fit1;
 end
%% subfunction for calculate fitness
 function [P]=deleteNeg(P)
%         P(:,index) = 0;
        [S,V]=find(abs(P)<=0.05);
        for i = 1:length(S)
            P(S(i),V(i))=0;
        end
        clear S V
        [S,V]=find(P>1); 
        for i = 1:length(S)
            P(S(i),V(i))=1;
        end
        [S,V]=find(P<-1); 
        for i = 1:length(S)
            P(S(i),V(i))=-1;
        end
 end
function re=notExist(ins, PosPop, Pop, NextPop, endl)
if ~isempty(PosPop)
for i=1:size(PosPop,1)              %in set PosPop
    if isequal(ins,PosPop(i,:)) == 1
        re = 0;
        return;
    end
end
end
for i=1:size(Pop,1)                %in set Pop
    if isequal(ins,Pop(i,:)) == 1
        re = 0;
        return;
    end
end
for i=endl:size(NextPop,1)                      %in set NextPop
    if isequal(ins,NextPop(i,:)) == 1
        re = 0;
        return;
    end
end
re = 1;
end


end