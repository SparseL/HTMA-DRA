
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
    lambda_max = norm( A'*b, 'inf' );
    lambda = 0.01*lambda_max;
    [w] = ISTA(A, b, lambda);
    
    networkLearn(:,nn) = w;
    fprintf('%dth node completed.\n', nn);
end
toc;
%%%%%%%  out of sample error  %%%%%%
%parameters
time = 10;
node = 20;         %node = n            
sample = 10;
data = zeros(time, node, sample);
datapredict = zeros(time, node, sample);
%data of t0
for spl = 1:sample
    for n = 1:node
        temp = rand;
        data(1,n,spl) = temp;
        datapredict(1,n,spl) = temp;
    end
end
%data based on network
for spl = 1:sample
    for ti = 1:(time-1)
        data(ti+1,:,spl) = sigmf(data(ti, :,spl)*network,[5 0]);
    end
end
%data based on networkPredict
for spl = 1:sample
    for ti = 1:(time-1)
        datapredict(ti+1,:,spl) = sigmf(data(ti,:,spl)*networkLearn,[5 0]);
    end
end
%out of sample error
res = abs(data - datapredict);
res = res.^2;
happen = sum(sum(sum(res)));
outofsampleerror = happen/(time-1)/node/sample
%fprintf('outofsampleerror is %f \n', outofsampleerror);
%%%%%%%  model error  %%%%%%%
modelerror = sum(sum(abs(networkLearn - network))) / (n*n)
%fprintf('modelerror is %f \n', modelerror);
%%%%%%%  ss mean  %%%%%%%%%%
network = abs(network);
networkPredict = abs(networkLearn);
FP = 0;
TP = 0;
TN = 0;
FN = 0;
for i = 1:n
    for j = 1:n
        if abs(networkPredict(i,j)) < 0.05
            if abs(network(i,j)) > 0.05
                FP = FP + 1;
            else 
                TP = TP + 1;
            end
        else
            if abs(network(i,j)) > 0.05
                TN = TN + 1;
            else
                FN = FN + 1;
            end
        end
    end
end
sensitivity = TP / (TP + FN)
specificity = TN / (TN + FP)
ss = (2 * specificity * sensitivity) / (specificity + sensitivity)

function [w] = ISTA(X, y, lambda, varargin)

maxIter = 100000;        %INPUT
optTol = 1e-5;
[n,p] = size(X);
t = (1-0.05)/norm(X'*X);
w = zeros(p,1);
k = p*0.2;                  %INPUT
for i = 1:maxIter
    [long,wide]=size(w);
    % Step in Negative Gradient Direction
    sig = X*w;
    w1 = w + t*X'*(y-sig);
    
    %lambdaNew
    sortB = sort(w1);
    KthB = abs(sortB(long+2-k));
    %lambdaNew=sqrt(96)/9*KthB^(3/2)*norm(X'*X);
    lambdaNew=sqrt(96)/9*((KthB^(3/2))/t);
    if lambda > lambdaNew
        lambda = lambdaNew;
    end
    
    % Soft Threshold
    condi= (54^(1/3)/4)*(lambda*t)^(2/3);
%     faiX = acos(lambda*t/8*(abs(w1/3).^(-3/2)));

    faiX = acos(lambda/8*(abs(w1/3).^(-3/2)));
    f = 2/3*w1.*(1+cos(2/3*pi-2/3*faiX));
    temp = max(abs(w1)-condi, 0);
    w2 = f.*sign(temp);   
    
    % Check for convergence
    if sum(abs(w-w2)) < optTol
        fprintf('Solution Found\n');
        break;
    end
    w = w2;
end
end
