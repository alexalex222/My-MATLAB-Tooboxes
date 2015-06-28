function [ W,P,Q,R,beta, XX_new, XY_new ] = RPLS2(XX_old, XY_old, x_new, y_new, A, lambda )
%Kernel algorithm of PLS
%   Based on Dayal & MacGregor 1997, Improved PLS algorithms
%	Author: Kuilin Chen
%	Date: 8-Jun-2013
%   Input
%   XX_old: old XX kernel matrix 
%   XY_old: old XY kernel matirx 
%   x_new: new x measurement
%   y_new: new y measurement
%   A: number of latent variables
%   
%   Output: 
%   W: weight matrix for defalted X
%   P: loading matrix for X
%   Q: loading matrix for Y
%   R: weight matrix for X
%   beta: regression ceofficients between X and Y
%   XX_new: new XX kernel matrix
%   XY_new: new XY kernel matrix

XX_new = lambda*XX_old + x_new'*x_new;
XY_new = lambda*XY_old + x_new'*y_new;

XX = XX_new;
XY = XY_new;

M = size(y_new,2);

% XY=X'*Y; % compute the covariance
% XX=X'*X; % matrices
for i=1:A, % A=number of PLS components to be computed
    if M == 1, % if there is a single response variable, compute the
    w = XY; % X-weights as shown here
    else % else
    [C,D]=eig(XY'*XY); % ?rst compute the eigenvectors of YTXXTX
    q=C(:,find(diag(D)==max(diag(D)))); % ?nd the eigenvector corresponding to the largest eigenvalue
    w=(XY*q); % compute X-weights
    end
    w=w/sqrt(w'*w); % normalize w to unity
    r=w; % loop to compute ri
    for j=1:i-1,
        r=r-(P(:,j)'*w)*R(:,j);
    end
    tt=(r'*XX*r); % compute tTt
    p=(r'*XX)'/tt; % X-loadings
    q=(r'*XY)'/tt; % Y-loadings
    XY=XY-(p*q')*tt; % XTY de?ation
    if i == 1
        W = w;
        P = p;
        Q = q;
        R = r;
    else
        W=[W w]; % storing loadings and weights
        P=[P p];
        Q=[Q q];
        R=[R r];
    end
end
    beta=R*Q';

end