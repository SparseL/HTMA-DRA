function [A,y,W] = data_generateFCM(SV,N,T,K)
A = zeros(K*SV,T);W = zeros(T,T);
W = sprandn(T,T,N);
for i = 1:T
    for j = 1:T
        if abs(W(i,j))>1
            W(i,j) = 2*rand-1;
        end
        if abs(W(i,j))>0 && abs(W(i,j))<0.05
            W(i,j) = 2*rand-1;
        end
    end
end

A = [];
for j = 1:SV
    a(1,:) = rand(1,T);
    for i = 2:K
        a(i, :) = sigmf(a(i-1,:)*W, [5 0]);
    end
    A = [A;a];
end
Aa = A;
for i = 1:SV
    A(i*K-(i-1),:) = [];
end
for i = 1:SV
    Aa((i-1)*K-(i-1)+1,:) = [];
end
y = -log((1-Aa)./Aa)/5;
% [s,v] = find(y == 0);
% for i = 1:length(s)
%     y(s(i),v(i)) = 0.000001;
% end
% [s1,v1] = find(y == 1);
% for i = 1:length(s1)
%     y(s1(i),v1(i)) = 0.999999;
% end

end