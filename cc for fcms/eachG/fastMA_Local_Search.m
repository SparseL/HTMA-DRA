%% editer: Kai Wu
%  opitimizaton: min ||x||_1 && ||AX-y||_2^2
%  x is a sparse singal. y is the observation. A is observation matrix.
%% main function 
function [P_best,P,plot_P_best,fit,spendbudget,fit2]...
    =fastMA_Local_Search(A,y,norm_index,MAX_interation,popsize) 
%% initialization
%-------------------initialization start-----------------------------------
rand('twister',sum(100*clock));
[~,N]=size(A);
spendbudget = 0;gen = 1;
% %initialize the population using lasso
% lambdaValue = [];P = [];
% lambdaValue = 10e.^(-6:2:-1); % the parameters need to be changed
lambdaValue = 10e-4;
% lambdaValue = [lambdaValue,1:10:100];
PP = lasso(A,y,'Lambda',lambdaValue,'Alpha',1,'MaxIter',10);
% PP = lasso(A,y,'Lambda',lambdaValue,'Alpha',1);
P = PP';

%-------------------initialization end--------------------------
%%%-----------------main loop start-----------------------------
while (gen<=MAX_interation)
% fprintf('CCMA generations: %d\n',gen);  
   
    %% local search   
%-----------local search----------------------------------------
    newP = [];
    for i = 1:popsize
        [newPP,lambda] = halflocalsearch(A, y, P(i,:));
        newP = [newP;newPP];
    end
%     newP=deleteNeg(newP); ii = 1;newi = [];
%     for i = 1:size(newP,1)-1
%         if notExist(newP(i,:),P,newP,i+1)==1        %avoid repetition
%             newi(ii) = i;
%             ii = ii+1;
%         end
%     end
%     newP_final = newP(newi,:);
    newP = [newP;P];
    [newPF,newPF2] = fitness(newP,lambda);
%     ii = 1;newi = [];
%     for j=1:size(newP_final,1)-1
%     for i=j+1:size(newP_final,1)                %in set Pop
%         if (abs(newP_fit(j,:)-newP_fit(i,:)) <= 10e-6)
%             re = 0;
%             return;
%         end
%     end
%     if re == 0
%         newi(ii) = j;
%         ii = ii+1;
%     end
%     end
%     newP = newP_final(newi,:);
%     newPF = newP_fit(newi);
%     newPF2 = newP_fit2(newi);
%     spendbudget = spendbudget + size(newP_final,1); 
%-----------------update population------------------------------------
    [~,LC]=sort(newPF,'ascend');
    P = newP(LC(1:popsize),:);
    fit = newPF(LC(1:popsize));
    fit2 = newPF2(LC(1:popsize));
%-----------------update best indidual-------------------------------------
    min_ = find(fit == min(fit));
    P_best = P(min_(1),:);
    plot_P_best(gen) = fit(min_(1));
    gen=gen+1;
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
        [S,V]=find(P<-1);
        for i = 1:length(S)
            P(S(i),V(i))=-1;
        end
        clear S V
        [S,V]=find(P>1); % the parameter needs to be optimized. 
        for i = 1:length(S)
            P(S(i),V(i))=1;
        end
 end
function re=notExist(ins, Pop, NextPop, endl)
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