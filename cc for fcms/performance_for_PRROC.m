function performance_for_PRROC
clc,clear
% 处理数据
Sv = 2:1:8;N = 100;pc=0.03;K=12;
fid_fg = fopen('100-k6-fg.txt','w');
fid_bg = fopen('100-k6-bg.txt','w');
for modelS = 1:3
for sizes = 1:4
tic
for data = 1:7
SV = Sv(data);
savefile=sprintf('modelS_%d_size_%d_N_%d_density_%d_SV_%d_learnt',modelS,sizes,N,pc*100,SV);
load(savefile,'-mat');
[S1,V1] = find(W>0);
[S2,V2] = find(W==0);

for kk = 1:5
xp_MA1 = xp_MA(:,:,kk);
xp_MA1 = topkprediction(xp_MA1,N,K);                                                                                                                                                                                                                                                                                                                   
for i = 1:length(S1)
    fg(1,i) = xp_MA1(S1(i),V1(i));
end
for i = 1:length(S2)
    bg(1,i) = xp_MA1(S2(i),V2(i));
end
fprintf(fid_fg,'%d\t',fg);
fprintf(fid_fg,'\n');
fprintf(fid_bg,'%d\t',bg);
fprintf(fid_bg,'\n');
clear fg bg
end
end
toc
end
end
fclose(fid_fg);
fclose(fid_bg);
end

function [Y] = topkprediction(X,N,K)
Y = zeros(N,N);
% [~,V] = sort(X,'descend'); 
% for i = 1:N
%     for j = 1:K
%         Y(V(j,i),i) = 1;
%     end
% end
[S11,V11] = find(X>=0.2);
[S12,V12] = find(X<0.2);
for i = 1:length(S11)
    Y(S11(i),V11(i)) = 1;
end
for i = 1:length(S12)
    Y(S12(i),V12(i)) = 0;
end
end