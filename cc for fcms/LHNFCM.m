% LHN-FCM
clc,clear
%No. of corresponding state vector squence
SV = 16;K = 6;
% W = load('supervisoryControl1.txt');
% W = load('Brazilian Amzaon Matrix.txt');
% W = load('EducationSoftwareMatrix.txt');
W = load('MPS.txt');
% n_zeros = [3,10,13,19,20,21,23];
T = size(W, 1);
% generate density N, nodes = T, datalength = SV*(K).
A = zeros(K*SV,T);
for j = 1:SV
    A((j-1)*K+1,:) = rand(1,T);
    for i = 2:K
        A((j-1)*K+i,:) = 1./(1+exp(-5*A((j-1)*K+i-1,:)*W));
    end
end

for j = 1:SV
    for i = 1:K-1
        y((j-1)*(K-1)+i,:) = -log((1-A((j-1)*(K)+i+1,:))./A((j-1)*(K)+i+1,:))/5;
    end
end

for i = 1:SV
    A(i*K-(i-1),:) = [];
end

Hnode = 3;
for j = 1:T
%     Wc(1,j) = find01(A, Hnode, j, W);
    if Hnode==j
        Wc(1,j) = inf;
        Wc1(1,j) = inf;
    else
        [Wc(1,j),Wc1(1,j)] = find01_LASSO(A, Hnode, j, W);
%         Wc(1,j) = find01(A, Hnode, j, W);
    end
end
figure;plot(Wc,'b-o');
xlabel('Node no.');
ylabel('\sigma');
set(gca,'FontSize',14);
set(get(gca,'YLabel'),'Fontsize',20);
set(get(gca,'XLabel'),'Fontsize',20);
