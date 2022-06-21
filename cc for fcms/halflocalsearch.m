function [solution,lambda] = halflocalsearch(X, y, w)
% case one: random selection from the input
% case two: the best k population
% the selected w, norm 0.
Kmin = 0; Kmax = 0.8*length(w);
e = zeros(size(w));
e(w~=0) = 1;
k1 = sum(e,2);
% T-neighborhood of k1.
T = 5;

if k1+T < Kmax
    kr = T+k1;
else
    kr = Kmax;
end

if k1-T > Kmin
    kL = k1-T;
else
    kL = Kmin;
end

if kr>kL
    k2 = kL+randperm(kr-kL,1);
else
    k2 = kL;
end

karchive = kL:1:k2;
% obtain Truncation solutions by multi-sparsity
[solution_k2,lambda] = halfL(X, y, w', k2);
solution = solution_k2';
[~,I] = sort(abs(solution_k2),'ascend');
for i = 1:length(karchive)-1
    temp = solution_k2';
    index = k2-karchive(i);
    temp(I(1:index)) = 0;
    solution = [solution;temp];
end
end

