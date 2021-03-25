clear all
clc
close all
cd('C:\Users\Zhang Yongfeng\Desktop\ZernikeCoefficient transformation\derivation')
%% Demonstration
syms a1 a2 b1 b2 real
A=[a1;a2];
B=[b1;b2];
direct=[];
indirect=[];
syms r real
for k=0:2:20
    direct=[direct;expand(DOT(A,B)^k)];
    SUM=sym(0);
    for r=0:floor(k/2)
        SUM=SUM+factorial(k)/factorial(r)/factorial(k-r)/...
            (1+(r==floor(k/2))*(mod(k,2)==0))*DOT(A,A)^r*DOT(B,B)^r*DOT(BINGM(A,k-2*r),BINGM(B,k-2*r))/2^(k-1);
    end
    indirect=[indirect;SUM];
end
diff=simplify(expand(direct-indirect))
