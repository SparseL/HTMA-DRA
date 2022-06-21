
% function FCM_L12
clear;
close all;
clc;

SV = 20;
% density
N = 0.2;
%Number of nodes
T = 100;
% number of data length
K = 10;
A = zeros(K*SV,T);W = zeros(T,T);

%Generate problem data
% network = load('EducationSoftwareMatrix.txt');
m = 10;       % number of examples                                        %INPUT
n = T;       % number of features                                        %INPUT
%p = 4/n;      % sparsity density                                         %INPUT
network = sprand(n,n,0.2);
for i = 1:n
    for j = 1:n
        if abs(network(i,j))>1
            network(i,j) = 2*rand-1;
        end
        if abs(network(i,j))>0 && abs(network(i,j))<0.05
            network(i,j) = 2*rand-1;
        end
    end
end
W = network;
for j = 1:SV
    A((j-1)*K+1,:) = rand(1,T);
    for i = 2:K
        A((j-1)*K+i,:) = 1./(1+exp(-5*(A((j-1)*K+i-1,:)-1)*W));
    end
end
[s,v] = find(A == 0);
for i = 1:length(s)
    A(s(i),v(i)) = 0.000001;
end
[s1,v1] = find(A == 1);
for i = 1:length(s1)
    A(s1(i),v1(i)) = 0.999999;
end
clear s v s1 v1
row = randperm(SV*K);

% observations
% y = A*W;
for j = 1:SV
    for i = 1:K-1
        y((j-1)*(K-1)+i,:) = -log((1-A((j-1)*(K)+i+1,:))./A((j-1)*(K)+i+1,:))/5;
    end
end

for i = 1:SV
    A(i*K-(i-1),:) = [];
end

tic;
for nn = 1:1
%     x0 = network(:,nn);
%     b = sigmf(A*x0, [5 0]) + sqrt(0.001)*rand(m,1);
%     b = A*x0 + sqrt(0.001)*rand(m,1);
    b = y(:,nn);
    lambda = 0.01;
    networkLearn(:,nn) = ISTA(A, b, lambda);  
    fprintf('%dth node completed.\n', nn);
end
toc;

%%%%%%%  ss mean  %%%%%%%%%%
% network = abs(network);
% networkPredict = abs(networkLearn);
% FP = 0;
% TP = 0;
% TN = 0;
% FN = 0;
% for i = 1:n
%     for j = 1:n
%         if abs(networkPredict(i,j)) < 0.05
%             if abs(network(i,j)) > 0.05
%                 FP = FP + 1;
%             else 
%                 TP = TP + 1;
%             end
%         else
%             if abs(network(i,j)) > 0.05
%                 TN = TN + 1;
%             else
%                 FN = FN + 1;
%             end
%         end
%     end
% end
% sensitivity = TP / (TP + FN)
% specificity = TN / (TN + FP)
% ss = (2 * specificity * sensitivity) / (specificity + sensitivity)

function [w] = ISTA(X, y, lambda, varargin)

maxIter = 1e6;        %INPUT
optTol = 1e-5;
[n,p] = size(X);
t = 0.5*(norm(X'*X)).^(-2);
w = zeros(p,1);
k = p*0.2;                  %INPUT
for i = 1:maxIter
    % Step in Negative Gradient Direction
    sig = X*w;
    w1 = w - t*X'*(sig-y);
    
    % Soft Threshold
    condi= (54^(1/3)/4)*(lambda*t)^(2/3);
    faiX = acos(lambda*t/8*(abs(w1/3).^(-3/2)));
    f = 2/3*w1.*(1+cos(2/3*pi-2/3*faiX));
    temp = max(abs(w1)-condi, 0);
    w2 = f.*sign(temp);   
    
    % Check for convergence
    if sum(abs(w-w2)) < optTol
        fprintf('Solution Found\n');
        break;
    end
    w = w2;
    [fit_(i),fit(i)]=fitness(w',y,X,lambda);
end
plot(fit_,'b-');
figure;plot(fit,'r-');
end

function [fit_,fit2]=fitness(P_,y,A,lambda)
        norm_index=1/2;
        [popsize_,~]=size(P_);
        fit1=zeros(popsize_,1);
        if norm_index==1
            fit1=sum(abs(P_));
        elseif norm_index==0
            e=zeros(size(P_));
            e(P_~=0)=1;
            fit1=sum(e,2);
        elseif norm_index==1/2
            temp=sum(sqrt(abs(P_)));
            fit1=temp.^2;
        end
        fit2=sum((abs(A*P_'-repmat(y,[1,popsize_]))).^2);
        % /(2*length(y))
         fit2= fit2';
         fit_ = fit1;
%          fit_ = fit2+lambda*fit1;
 end
