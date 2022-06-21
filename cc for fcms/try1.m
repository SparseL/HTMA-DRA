clc,clear
nf = 0.0;
%No. of corresponding state vector squence
SV = 2;
% density
N = 0.2;
%Number of nodes
T = 20;
% number of data length
K = 10;
% generate density N, nodes = T, datalength = SV*(K).
A = zeros(K*SV,T);
W = zeros(T,T);
n = 0;

% W = load('supervisoryControl.txt');
% T = size(W, 1);
% W = -1+2*rand(T,T);
for realization = 1:1
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

for j = 1:SV
    A((j-1)*K+1,:) = rand(1,T);
    for i = 2:K
        A((j-1)*K+i,:) = 1./(1+exp(-5*(A((j-1)*K+i-1,:)-1)*W));
    end
end
[s,v] = find(A == 0);
for i = 1:length(s)
    A(s(i),v(i)) = 0.000001;
end
[s1,v1] = find(A == 1);
for i = 1:length(s1)
    A(s1(i),v1(i)) = 0.999999;
end
clear s v s1 v1
row = randperm(SV*K);

% observations
% y = A*W;
for j = 1:SV
    for i = 1:K-1
        y((j-1)*(K-1)+i,:) = -log((1-A((j-1)*(K)+i+1,:))./A((j-1)*(K)+i+1,:))/5;
    end
end

for i = 1:SV
    A(i*K-(i-1),:) = [];
end

tic
for i = 1:T
    x0 = A\y(:,i);
    xp(:,i) = l1eq_pd(x0, A, [], y(:,i), 1e-3);
end
toc

for i = 1:T
    for j = 1:T
        if abs(xp(i,j)) <= 0.05
             xp(i,j) = 0;
        end
     end
end

tp = 0;
tn = 0;
fn = 0;
fp = 0;
% calculate model error
for i = 1:T
    for j = 1:T
        if abs(xp(i,j)) >= 0.05 && abs(W(i,j)) >= 0.05
            tp = tp+1;
        elseif abs(xp(i,j)) < 0.05 && abs(W(i,j)) < 0.05
            tn = tn+1;
        elseif abs(xp(i,j)) >= 0.05 && abs(W(i,j)) < 0.05
            fn = fn+1;
        else
            fp = fp+1;
        end
    end
end
specificity(realization) = tp/(tp+fn);
sensitivity(realization) = tn/(tn+fp);
SS_mean(realization) = 2*specificity(realization)*sensitivity(realization)/(sensitivity(realization)+specificity(realization));
disp(sprintf('SS mean = %8.3e,specificity = %8.3e,sensitivity = %8.3e', SS_mean,specificity,sensitivity));
%fprintf(fp1,'%f\t%f\t%f\t',SS_mean,specificity,sensitivity);

% % compute data error
%  for j = 1:SV
%      AA((j-1)*K+1,:) = A((j-1)*K+1,:);
%     for i = 2:K
%         AA((j-1)*K+i,:) = 1./(1+exp(-5*A((j-1)*K+i-1,:)*xp));
%     end
%  end
% yy = AA*xp;
% data_error(realization) = sum(sum((y-yy).*(y-yy)))/(SV*T*K);
% 
% % compute model error
% model_error(realization) = sum(sum(abs(xp-W)))/(T*T);
% 
% % compute out of sample error
% for i = 1:10
%     for j = 1:SV
%         out_error((j-1)*K+1,:) = rand(1,T);
%         for i = 2:K
%              out_error((j-1)*K+i,:) = 1./(1+exp(-5* out_error((j-1)*K+i-1,:)*W));
%         end
%     end
%     Y = out_error*W;
%    for j = 1:SV
%         Aa((j-1)*K+1,:) = out_error((j-1)*K+1,:);
%         for i = 2:K
%             Aa((j-1)*K+i,:) = 1./(1+exp(-5* out_error((j-1)*K+i-1,:)*xp));
%         end
%     end
%     YY = Aa*xp;
%     out_of_sample_error(i,:) = sum(sum(abs(Y-YY)))/(T*SV*K);
% end
% out1(realization) = mean(out_of_sample_error);
% end
% disp(sprintf('SS mean = %8.3e,specificity = %8.3e,sensitivity = %8.3e', mean(SS_mean),mean(specificity),mean(sensitivity)));
% disp(sprintf('data error = %8.3e, model error = %8.3e', mean(data_error),mean(model_error)));
% disp(sprintf('data error_std = %8.3e, model error_std = %8.3e', std(data_error),std(model_error)));
% disp(sprintf('out_mean = %8.3e, out_std = %8.3e', mean(out1),std(out1)));
% % row = randi(T, 4, 1);
% row = randperm(T);
% tic
% for i = 1:10
%     for j = 1:T
%         Wc(i,j) = find01(A, row(i), j, W);
%     end
% end
% toc
% % for i = 1:10
% %      for j = 1:T
% %           if abs(Wc(i,j)) >= 0.05 && abs(W(row(i),j)) >= 0.05
% %               tn = tn+1;
% %           elseif abs(Wc(i,j)) < 0.05 && abs(W(row(i),j)) < 0.05
% %               tp = tp+1;
% %           elseif abs(Wc(i,j)) >= 0.05 && abs(W(row(i),j)) < 0.05
% %               fn = fn+1;
% %           else
% %               fp = fp+1;
% %           end    
% %      end
% % end

tic
for i = 1:T
    for j = 1:T
        if i ~= j
            [Wc(i,j), Wc1(i,j)] = find01(y, A, i, j);
        else
            Wc(i,j) = 0;
            Wc1(i,j) = 0;
        end
        
    end
end
toc

% for i = 1:T
%     for j = 1:T
%         if abs(Wc(i,j)) >= 0.05 && abs(W(i,j)) >= 0.05
%            tn = tn+1;
%         elseif abs(Wc(i,j)) < 0.05 && abs(W(i,j)) < 0.05
%             tp = tp+1;
%         elseif abs(Wc(i,j)) >= 0.05 && abs(W(i,j)) < 0.05
%             fn = fn+1;
%         else
%             fp = fp+1;
%         end
%     end
% end
% % specificity = tp/(tp+fn)
% % sensitivity = tn/(tn+fp)
% % SS_mean = 2*specificity*sensitivity/(sensitivity+specificity)

% k = 1;
% kk = 1;
% for i = 1:T
%     for j = 1:T
%         if abs(W(i, j))>=0.05
%             existlink(1,k) = Wc(i, j);
%             k = k+1;
%         else
%             nulllink(1,kk) = Wc(i, j);
%             kk = kk+1;
%         end
%         
%     end
% end

k1 = 1;
for k = 0:0.01:1
    tp = 0;
    tn = 0;
    fn = 0;
    fp = 0;
    for i = 1:T
        for j = 1:T
            if i ~= j
            if abs(Wc(i,j)) >= k && abs(W(i,j)) >= 0.05
                tn = tn+1;
            elseif abs(Wc(i,j)) < k && abs(W(i,j)) < 0.05
                tp = tp+1;
            elseif abs(Wc(i,j)) >= k && abs(W(i,j)) < 0.05
                fn = fn+1;
            else
                fp = fp+1;
            end    
            end
        end
    end
    specificity(k1) = tp/(tp+fn);
    sensitivity(k1) = tn/(tn+fp);
    SS_mean(k1) = 2*specificity(k1)*sensitivity(k1)/(sensitivity(k1)+specificity(k1));
    k1 = k1+1;
end

% k1 = 1;
% for k = 0:0.01:1
%     for i = 1:10
%         for j = 1:T
%             if abs(Wc(i,j)) >= k && abs(W(row(i),j)) >= 0.05
%                 tn = tn+1;
%             elseif abs(Wc(i,j)) < k && abs(W(row(i),j)) < 0.05
%                 tp = tp+1;
%             elseif abs(Wc(i,j)) >= k && abs(W(row(i),j)) < 0.05
%                 fn = fn+1;
%             else
%                 fp = fp+1;
%             end    
%         end
%     end
%     specificity(k1) = tp/(tp+fn);
%     sensitivity(k1) = tn/(tn+fp);
%     SS_mean(k1) = 2*specificity(k1)*sensitivity(k1)/(sensitivity(k1)+specificity(k1));
%     k1 = k1+1;
% end

bestSSmean(realization) = max(SS_mean)
[c, v] = find(SS_mean == bestSSmean(realization));
sigma = v/100

% [c, v] = find(W ~= 0);
% for i = 1:length(c)
%      xx1(i) = W(c(i),v(i));
% end
% scatter(c,v,25,xx1,'full');
% colorbar
% xlabel('\itN');
% ylabel('\itN');
% figure
% [c1, v1] = find(abs(Wc) > 0.05);
% for i = 1:length(c1)
%      xx2(i) = Wc(c1(i),v1(i));
% end
% scatter(c1,v1,25,xx2,'full');
% colorbar
% xlabel('\itN');
% ylabel('\itN');
% for i = 1:T
%     for j = 1:T
%         if Wc(i,j)>sigma
%             Wc(i,j) = 1;
%         else
%             Wc(i, j) = 0;
%         end
%     end
% end
% [c2, v2] = find(abs(Wc) > 0.05);
% figure
% scatter(c2,v2,25,'full');
% xlabel('\itN');
% ylabel('\itN');

% hold on
% plot(SS_mean, 'g*-');
% axis tight
% xlabel('\sigma');
% ylabel('SS mean');
% set(gca,'xticklabel',[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]);
end