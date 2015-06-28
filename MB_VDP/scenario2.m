clear;
clc;
close all;
load TEPdata\mode1_normal;
load TEPdata\mode1_idv1;
load TEPdata\mode1_idv2;
load TEPdata\mode1_idv3;
load TEPdata\mode1_idv4;
load TEPdata\mode1_idv5;
load TEPdata\mode1_idv6;
load TEPdata\mode1_idv7;
load TEPdata\mode1_idv8;
load TEPdata\mode1_idv9;
load TEPdata\mode1_idv10;
load TEPdata\mode1_idv11;
load TEPdata\mode1_idv12;
load TEPdata\mode1_idv13;
load TEPdata\mode1_idv14;
load TEPdata\mode1_idv15;

load TEPdata\mode3_normal;
load TEPdata\mode3_idv1;
load TEPdata\mode3_idv2;
load TEPdata\mode3_idv3;
load TEPdata\mode3_idv4;
load TEPdata\mode3_idv5;
load TEPdata\mode3_idv6;
load TEPdata\mode3_idv7;
load TEPdata\mode3_idv8;
load TEPdata\mode3_idv9;
load TEPdata\mode3_idv10;
load TEPdata\mode3_idv11;
load TEPdata\mode3_idv12;
load TEPdata\mode3_idv13;
load TEPdata\mode3_idv14;
load TEPdata\mode3_idv15;

mode1n=[mode1_normal.x,mode1_normal.y];
mode3n=[mode3_normal.x,mode3_normal.y];
mode1f=[mode1_idv8.x,mode1_idv8.y];
mode3f=[mode3_idv13.x,mode3_idv13.y];
% mode1f=[mode1_idv2.x,mode1_idv2.y;mode1_idv12.x,mode1_idv12.y;mode1_idv13.x,mode1_idv13.y;mode1_idv14.x,mode1_idv14.y;];
% mode1f=mode1f([1:4:4000],:);
Tdata=[mode1n;mode3n;mode3f];


[nr,nc]=size(mode1n);
[nr1,nc1]=size(Tdata);
[~,mu,sigma]=zscore(mode1n);
mu=repmat(mu,nr1,1);
sigma=repmat(sigma,nr1,1);
dataP1=(Tdata-mu)./sigma;

[T1,P1] = pca(dataP1,10);

% X=[mode1_normal.x;mode3_normal.x;mode1_idv1.x];
% Y=[mode1_normal.y;mode3_normal.y;mode1_idv1.y];
% Xz=zscore(X);
% Yz=zscore(Y);
% [~,~,Xs,Ys]=plsregress(Xz,Yz,7);
% T1=[Xs,Ys];



results = MB_VDP(T1',3000,0,generate_prior(T1'));
[~,classlabel]=max(results.q_z.singlets,[],2);
for i=1:1000
    classlabel(i)=1;
end

for i=1001:3000
    classlabel(i)=2;
end
figure;
scatterMixture(T1(:,1:2),classlabel);
legend('Base Mode','New Mode','Outlier','Outlier');
xlabel('PC1');
ylabel('PC2');
% legend('Base Mode','New Mode');