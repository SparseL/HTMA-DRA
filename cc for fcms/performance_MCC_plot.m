% performance_MCC
function performance_MCC_plot
clc,clear
% 处理数据
Sv = 2:1:8;N=100;pc=0.03;
fid_MCC = fopen('k3mcc.txt','w');
tic
for modelS = 1:3
for sizes = 1:4
for data = 1:7
SV = Sv(data);

% savefile=sprintf('N_%d_SV_%d_modelS_%d_lambda_%d_learnt',N,SV,modelS,la_step);
% load(savefile,'-mat');

MCC = zeros(1,5);
savefile=sprintf('modelS_%d_size_%d_N_%d_density_%d_SV_%d_learnt',modelS,sizes,N,pc*100,SV);
load(savefile,'-mat');

for kk = 1:5
xp_MA2 = xp_MA(:,:,kk);
xp_MA2 = topkprediction(xp_MA2,N);                                                                                                                                                                                                                                                                                                                   
MCC(1,kk) = MCC_predict(W,xp_MA2);

end
fprintf(fid_MCC,'%d\t',MCC');
fprintf(fid_MCC,'\n');
end
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