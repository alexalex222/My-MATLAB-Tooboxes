clear java;
clear;
close all;

javaaddpath(pwd);

X = rand(100,5);
y = rand(100,1);
X1  = [ones(100,1) X];
mdl = LinearModel.fit(X,y)


regression = MultipleLinearRegression(X,y,true);

beta = regression.getbeta();
betaSD = regression.getbetaSD();
tstats = regression.gettStats();
R2 = regression.getR2();
R2bnar = regression.getR2bar;
s2 = regression.gets2();
Fstats = regression.getFstats();


