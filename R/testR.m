clear;
close all;
clc;

% openR; %Open connection to R server 
% x=[1:50]; %create x values in Matlab 
% putRdata('x = as.vector(x)',x); %put data into R workspace 
% evalR('y<-sqrt(x)'); %evaluate in R 
% evalR('plot(x,y)') %plot in R 

load hald
X = ingredients; % Predictor variables
y = heat; % Response

% mdl = LinearModel.fit(X,y)

% openR;
putRdata('x1', ingredients(:,1));
putRdata('x2', ingredients(:,2));
putRdata('x3', ingredients(:,3));
putRdata('x4', ingredients(:,4));
evalR('x1 = as.vector(x1)');
evalR('x2 = as.vector(x2)');
evalR('x3 = as.vector(x3)');
evalR('x4 = as.vector(x4)');

putRdata('y', heat);
evalR('y = as.vector(y)');
% evalR('library(car)');
% evalR('head(Moore)');
% [result,status,msg]=evalR('model <- lm(formula = conformity ~ partner.status*fcategory, data = Moore)');
[result,status,msg] = evalR('model <- lm(y ~ x1 + x2 + x3 + x4)',0);
[c,status,msg] = evalR('c <- model$coef[[1]]',0); error(msg)
[m,status,msg] = evalR('m <- model$coef[[2]]',0); error(msg)
evalR('layout(matrix(c(1,2,3,4),2,2))');
evalR('plot(model)');
[interceptR,status,msg] = getRdata('c'); error(msg)
[slopeR,status,msg] = getRdata('m'); error(msg)
% [Rlm_model, status, msg] = getRdata('ft1');