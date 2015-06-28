function pls_cv_100out = plscv100out(x,y,vl)
% 
% PLS cross validation leave-100-out
%

[sample,dim]=size(x);
folder=fix(sample/100);

RMSE=zeros(folder,1);

for i=1:folder
    e=i*100;
    s=e-99;
    % Define crossvalidation part
    x_val=x(s:e,:);
    y_val=y(s:e,:);
    
    % Define training part
    x_train=x([1:s-1,e+1:end],:);
    y_train=y([1:s-1,e+1:end],:);
    
    [XL,YL,XS,YS,BETA] = plsregress(x_train,y_train,vl);
    
    y_hat = [ones(100,1),x_val]*BETA;
    
    RMSE(i)= mean((y_hat-y_val).^2);
end

pls_cv_100out=mean(RMSE);

    
    
