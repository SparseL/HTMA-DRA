clc,clear 
%No. of corresponding state vector squence
SV = 10;
% density
N = 0.2;
%Number of nodes
T = 40;
% number of data length
K = 6;
% generate density N, nodes = T, datalength = SV*(K).
[A,y,W] = data_generateFCM(SV,N,T,K);