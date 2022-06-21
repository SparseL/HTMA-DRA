n = 512;%n is the dimension of the signal x
d=330;
sparsityRatio=0.05;
noiseVariance=0;
% k=ceil(sparsityRatio*n);
k = 30;
%signal=1/randn
signalname='randn'; %||A*x-y||_2^2
sparseSupport = randperm(n); x0=zeros(n,1);
x0(sparseSupport(1:k))=randn(1,k);
normx0=norm(x0,2); x0 = x0 / norm(x0);
k=length(nonzeros(x0));
% Generate Gaussian dictionary
AMatrix = randn(d,n); matrixNorm = AMatrix.'*AMatrix;
matrixNorm = sqrt(diag(matrixNorm)).';
AMatrix = AMatrix./repmat(matrixNorm, [d,1]);
noise_added=noiseVariance*randn(d,1);
y = AMatrix*x0; y = y +noise_added;
[w,gen] = halfL(AMatrix, y, [], k);
MSE = norm(w-x0);