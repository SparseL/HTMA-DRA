% for i = 1:40
%     min_ = find(fit(:,i) == min(fit(:,i)));
%     plot_P(i) = fit(min_(1));
%     min1 = find(fitini(:,i) == min(fitini(:,i)));
%     plot_P1(i) = fit(min1(1));
% end
% figure
% plot(plot_P1-plot_P,'b-o')
% 
% xlabel('No. Node','FontName','Times new roman','fontsize',12,'FontWeight','bold');
% ylabel('Fitness \itf','FontName','Times new roman','fontsize',12,'FontWeight','bold');
% legend({'With DRAS','Without DRAS'},'FontName','Times new roman','FontSize',12,'FontWeight','bold');
% clc,clear
SV = 20;N = 0.2;T = 40;K = 6;
savefile=sprintf('data-Node_%d_density_%d_SV_%d',T,N*100,SV);
load(savefile);
[AA,yy] = dataG(W,10,K,T);
for tryindex = 1:1

    savefile=sprintf('CCMA-Node_%d_density_%d_SV_%d',T,N*100,SV);
    load(savefile);
    savefile=sprintf('CommonCC-Node_%d_density_%d_SV_%d',T,N*100,SV);
    load(savefile);
    [data(tryindex),out(tryindex),model(tryindex),SS(tryindex)] = measureFCM(A,y,xp_MA,SV,T,K,W,AA,yy);
    [data1(tryindex),out1(tryindex),model1(tryindex),SS1(tryindex)] = measureFCM(A,y,xp_MA1,SV,T,K,W,AA,yy);
end

function [A,y] = dataG(W,SV,K,T)
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
end