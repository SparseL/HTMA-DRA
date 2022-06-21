function [P,fit,fit2]=MA_Initialization(A,y,norm_index,popsize,lambda) 
%% initialization
%-------------------initialization start-----------------------------------

[~,N]=size(A);
P=zeros(popsize,N);
% %initialize the population using lasso
% lambdaValue = [];P = [];
lambdaValue = 10.^(-6:2:-1); % the parameters need to be changed
% lambdaValue = [lambdaValue,1:10:100];
PP = lasso(A,y,'Lambda',lambdaValue,'Alpha',1,'MaxIter',30);
% PP = lasso(A,y,'Lambda',lambdaValue,'Alpha',1);
P = PP';
remainP = popsize-length(lambdaValue);
W = -1+2*rand(remainP,N);

% W = sprandn(popsize,N,0.3);
for i1 = 1:remainP
    for j1 = 1:N
        if rand < 0.3
            W(i1,j1) = 0;
        end
    end
end
P = [P;W];
[P]=deleteNeg(P);
% P=deleteNeg(P,index);
[fit,fit2]=fitness(P,lambda);


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
        fit2=sum((abs(A*P_'-repmat(y,[1,popsize_]))).^2,1);
        % /(2*length(y))
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
end
