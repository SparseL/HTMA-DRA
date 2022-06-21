% performance_MCC
function performance_MCC
clc,clear
% 处理数据
Sv = 2:1:8;N=100;pc=0.06;
fid_MCC = fopen('lambda-CN-MCC.txt','w');
lambdaV = [10,1,0.1,0.01,0.001,0.0001,0.00001];
tic
for sizes = 1:4
for data = 3:3
SV = Sv(data);
MCC = zeros(7,5);
MCC_Mean = zeros(7,2);
% savefile=sprintf('N_%d_SV_%d_modelS_%d_lambda_%d_learnt',N,SV,modelS,la_step);
% load(savefile,'-mat');
for la_step = 1:length(lambdaV)
savefile=sprintf('size_%d_N_%d_density_%d_SV_%d_lambda_%d_learnt',sizes,N,pc*100,SV,la_step);
load(savefile,'-mat');

for kk = 1:5
xp_MA2 = xp_MA(:,:,kk);
xp_MA2 = topkprediction(xp_MA2,N);                                                                                                                                                                                                                                                                                                                   
MCC(la_step,kk) = MCC_predict(W,xp_MA2);
end

end
MCC_Mean(:,1) = mean(MCC,2);
MCC_Mean(:,2) = std(MCC');
fprintf(fid_MCC,'%d\t',MCC_Mean(:,1)');
fprintf(fid_MCC,'\n');
fprintf(fid_MCC,'%d\t',MCC_Mean(:,2)');
fprintf(fid_MCC,'\n');
end
end
fclose(fid_MCC);
end

function [Y] = topkprediction(X,N)
Y = zeros(N,N);
[S11,V11] = find(X>=0.2);
[S12,V12] = find(X<0.2);
for i = 1:length(S11)
    Y(S11(i),V11(i)) = 1;
end
for i = 1:length(S12)
    Y(S12(i),V12(i)) = 0;
end
end

function [MCC] = MCC_predict(x,y)
TP = 0;FP = 0; FN = 0;TN = 0;
le = size(x,1);
for i = 1:le
    for j = 1:le
        if x(i,j)>0 && y(i,j)==1
            TP = TP+1;
        elseif x(i,j)>0 && y(i,j)==0
            FN = FN+1;
        elseif x(i,j)==0 && y(i,j)==0
            TN = TN+1;
        else
            FP = FP+1;
        end
    end
end

MCC = (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
end