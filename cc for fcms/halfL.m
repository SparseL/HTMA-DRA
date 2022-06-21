function [w,gen] = halfL(X, y, w, k)
n = size(X,2);
if isempty(w)
    w = zeros(n,1);
end
maxIter = 20;        %INPUT
gen = zeros(maxIter,1);
optTol = 1e-8;
% t = (1-0.001)/norm(X'*X);
mtemp = sum(sum(X.^2));
t = (1-0.001)/mtemp;
 % input sparsity k norm(X'*X)
old_lambda = inf;
for i = 1:maxIter
%     gen(i) = sum((X*w-y).^2);
    % Step in Negative Gradient Direction
    w1 = w - t*X'*(X*w-y);
    
    %lambdaNew
    sortB = sort(w1,'descend');
    KthB = abs(sortB(k+1));
%     lambda=sqrt(96)/9*KthB^(3/2)*norm(X'*X);
    lambda=sqrt(96)/9*KthB^(3/2)*mtemp;
    lambda = min(old_lambda,lambda);
    old_lambda = lambda;
    % half Threshold
    condi= ((54^(1/3))/4)*(lambda*t)^(2/3);
    faiX = acos(lambda*t/8*(abs(w1/3).^(-3/2)));
    f = (2/3)*w1.*(1+cos((2/3)*pi-(2/3)*faiX));
    temp = max(abs(w1)-condi, 0);
    w2 = f.*sign(temp);   
    
    % Check for convergence
    if sum(abs(w-w2)) < optTol
        fprintf('Solution Found\n');
        break;
    end
    w = w2;
end
