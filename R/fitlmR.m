function [ lm ] = fitlmR( X, y, intercept)
%multiple linear regression based on R interface
%   Input:
%   X: regressor matrix
%   y: response variable
%   intercept: flag to include intercept
%   Output:
%   lm: linear model structure

if nargin<2 || nargin>3
    error('2 or 3 inputs only')
end

if (nargin ==2)
    intercept = true;
end 


openR;

putRdata('X', X);
evalR('X = as.matrix(X)',0);
evalR('X = data.frame(X)',0);

putRdata('y', y);
evalR('y = as.vector(y)',0);
if (intercept)
    evalR('model <- lm(y ~ ., data = X)',0);
else
    evalR('model <- lm(y ~ .-1, data = X)',0);
end

evalR('summarylm <- summary(model)',0)
evalR('Fstats <- summarylm$fstatistic',0);
Fstats = getRdata('Fstats');

evalR('df <- summarylm$df',0);
df = getRdata('df');

evalR('resStd <- summarylm$sigma',0);
resStd = getRdata('resStd');

evalR('temptable <- summarylm$coefficients',0);
temptable = getRdata('temptable');
beta = temptable(:,1);
betaStdErr = temptable(:,2);
tStats = temptable(:,3);
pValue = temptable(:,4);

evalR('R2 <- summarylm$r.squared',0);
R2 = getRdata('R2');

evalR('R2A <- summarylm$adj.r.squared',0);
R2A = getRdata('R2A');

evalR('y_hat <- fitted(model)',0);
y_hat = getRdata('y_hat');
y_hat = y_hat';

evalR('y_res <- residuals(model)',0);
y_res = getRdata('y_res');
y_res = y_res';

evalR('CI <- confint(model, level=0.95)',0);
CI = getRdata('CI');

% plot regression results
% evalR('plot(y, y_hat, xlab ="y", ylab="y_hat")',0);
% evalR('abline(lm(y_hat~y))',0);

% plot regression results
evalR('plot(y, y_hat, xlab ="y", ylab="y_hat")',0);

% scatter plot
evalR('windows()',0);
evalR('y = data.frame(y)',0);
evalR('all <- cbind(X,y)',0);
evalR('pairs(all)',0);

% plot regression residuals
evalR('windows()',0);
evalR('layout(matrix(c(1,2,3,4),2,2))');
evalR('plot(model)');


lm.y_hat = y_hat;
lm.y_res = y_res;
lm.beta = beta;
lm.CI = CI;
lm.betaStdErr = betaStdErr;
lm.tStats = tStats;
lm.pValue = pValue;
lm.Fstats = Fstats;
lm.df = df;
lm.R2 = R2;
lm.R2A = R2A;
lm.resStd = resStd;


end

