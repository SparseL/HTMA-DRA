function [data_error,out_error,model_error,SS_mean] = measureFCM(A,y,xp,SV,T,K,W,AA,yy)
for realization = 1:1
for i = 1:T
    for j = 1:T
        if abs(xp(i,j)) <= 0.05
             xp(i,j) = 0;
        end
     end
end
% compute data error
data_error(realization) = sum(sum((abs(A*xp-y)).^2))/(SV*T*(K-1));
% compute out of sample error
% [AA,yy] = dataG(W,10,K,T);

out_error(realization) = sum(sum((abs(AA*xp-yy)).^2))/(10*T*(K-1));
% compute model error
model_error(realization) = sum(sum(abs(xp-W)))/(T*T);
for i = 1:T
    for j = 1:T
        if abs(xp(i,j)) <= 0.05
             xp(i,j) = 0;
        end
     end
end
tp = 0;tn = 0;
fn = 0;fp = 0;
% calculate SS mean
for i = 1:T
    for j = 1:T
        if abs(xp(i,j)) > 0 && abs(W(i,j)) > 0
           tn = tn+1;
        elseif abs(xp(i,j)) == 0 && abs(W(i,j)) == 0
            tp = tp+1;
        elseif abs(xp(i,j)) > 0 && abs(W(i,j)) == 0
            fn = fn+1;
        else
            fp = fp+1;
        end
    end
end
sensitivity(realization) = tn/(tn+fp);
specificity(realization) = tp/(tp+fn);
SS_mean(realization) = 2*specificity(realization)*sensitivity(realization)...
    /(sensitivity(realization)+specificity(realization));
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
