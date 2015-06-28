function [T,P] = pca(X,nPC)

% PCA: Principal Component Analysis
% NIPALS: Non-linear Iterative Projections using Alternating Least Squares

[m,n] = size(X);

% allocate T and P
T=zeros(m,nPC);
P=zeros(nPC,n);
% tol for convergence
tol=sqrt(eps);
for k=1:nPC
    %find the column which has the maximum norm
    [dum,idx]=max(sum(X.*X));
    t=X(:,idx);
    %storage to judge convergence
    t0=t-t;
    %iteration if not converged
    while norm(t-t0)>tol
        %iteration to approach the eigenvector direction
        p=X'*t;
        %normalize the vector
        p=p/norm(p);
        %save previous t
        t0=t;
        %t is a product of eigenvalue and eigenvector
        t=X*p;
    end
    %subtracing PC identified
    X=X-t*p';
    T(:,k)=t;
    P(k,:)=p;
end