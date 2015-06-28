clear;
clc;
close all;

load carbig;
X = [Displacement Horsepower Weight Acceleration MPG];
nans = sum(isnan(X),2) > 0;
X1=X(~nans,1:3);
Y1=X(~nans,4:5);
X1=zscore(X1);
Y1=zscore(Y1);
[A, B, r] = cca(X1,Y1);

plot(U(:,1),V(:,1),'.')
xlabel('0.0025*Disp+0.020*HP-0.000025*Wgt')
ylabel('-0.17*Accel-0.092*MPG')


% X1_cen=X1-repmat(mean(X1),size(X1,1),1);
% Y1_cen=Y1-repmat(mean(Y1),size(Y1,1),1);

I1=A'*(X1'*X1)*A;
I2=B'*(Y1'*Y1)*B;