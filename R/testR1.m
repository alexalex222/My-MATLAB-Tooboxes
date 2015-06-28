clear;
clc;
% Generate some synthetic linear regression data
n = 100; x = (1:n); y = 3*x + 500 + 30*randn(1,n);
% Do the linear regression in R
[status,msg] = putRdata('x',x); error(msg)
[status,msg] = putRdata('y',y); error(msg)
 [result,status,msg] = evalR('x = as.vector(x)',0);
  [result,status,msg] = evalR('y = as.vector(y)',0);
[result,status,msg] = evalR('model <- lm(y ~ x)',0); error(msg)
[result,status,msg] = evalR('c <- model$coef[[1]]',0); error(msg)
[result,status,msg] = evalR('m <- model$coef[[2]]',0); error(msg)
% Get the results back in Matlab
[interceptR,status,msg] = getRdata('c'); error(msg)
[slopeR,status,msg] = getRdata('m'); error(msg)