function plot_PRROC
clc,clear
% a = load('100-k6-AUROC.txt');
a = load('100-k6-AUPR.txt');

for sizes = 1:4
figure
% subplot(2,2,sizes);
% bar(a((sizes-1)*9+1:(sizes-1)*9+9,:))
% h=plot(a((sizes-1)*9+1,:),'*r-');set(h,'linewidth',1.5);
h=plot(a((sizes-1)*7+1:sizes*7,1),'*r--');set(h,'linewidth',1.5);
axis tight
hold on
h=plot(a((sizes-1)*7+1:sizes*7,2),'ob--');set(h,'linewidth',1.5);
h=plot(a((sizes-1)*7+1:sizes*7,3),'pg--');set(h,'linewidth',1.5);
h=plot(a((sizes-1)*7+1:sizes*7,4),'+y--');set(h,'linewidth',1.5);
h=plot(a((sizes-1)*7+1:sizes*7,5),'xm--');set(h,'linewidth',1.5);
h=plot(a((sizes-1)*7+1:sizes*7,6),'dc--');set(h,'linewidth',1.5);
h=plot(a((sizes-1)*7+1:sizes*7,7),'hk--');set(h,'linewidth',1.5);

h=plot(a((sizes-1)*7+1:sizes*7,8),'<c--');set(h,'linewidth',1.5);
h=plot(mean(a((sizes-1)*7+1:sizes*7,9:13),2),'sk-');set(h,'linewidth',1.5);

xlabel('N_M','FontName','Times new roman','fontsize',12,'FontWeight','bold');
ylabel('AUPR','FontName','Times new roman','fontsize',12,'FontWeight','bold');
legend({'BP','FIST','Homotopy','L1LS','LARS','LASSO','OMP','PALM','MASTNet'},'FontName','Times new roman','FontSize',12,'FontWeight','bold');
% legend({'N_M=2','N_M=3','N_M=4','N_M=5','N_M=6','N_M=7','N_M=8'},'FontName','Times new roman','FontSize',12,'FontWeight','bold')
% set(gca,'XTickLabel',{'BP','FIST','Homotopy','L1LS','LARS','LASSO','OMP','PALM','MASTNet'});
% legend({'BP','FIST','L1LS','LARS','LASSO','PALM','MASTNet'},'FontName','Times new roman','FontSize',12,'FontWeight','bold');
set(gca,'FontName','Times new roman','XTickLabel',{'N_M=2','N_M=3','N_M=4','N_M=5','N_M=6','N_M=7','N_M=8'});
end
% legend({'N_M=2','N_M=3','N_M=4','N_M=5','N_M=6','N_M=9','N_M=8'},'FontSize',12,'FontWeight','bold')
