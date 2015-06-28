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

% load bcmode1idv8;
% classlabel1=classlabel;
load bcmode3mix;
% classlabel2=classlabel;

mode1n=[mode1_normal.x,mode1_normal.y];
mode3n=[mode3_normal.x,mode3_normal.y];
mode1f=[mode1_idv8.x,mode1_idv8.y];
% mode3f=[mode3_idv13.x,mode3_idv13.y];
% mode1f=[mode1_idv2.x,mode1_idv2.y;mode1_idv12.x,mode1_idv12.y;mode1_idv13.x,mode1_idv13.y;mode1_idv14.x,mode1_idv14.y;];
% mode1f=mode1f([1:4:4000],:);
mode3f=[mode3_idv1.x,mode3_idv1.y;mode3_idv10.x,mode3_idv10.y;mode3_idv13.x,mode3_idv13.y;mode3_idv14.x,mode3_idv14.y;];
mode3f=mode3f([1:4:4000],:);
Tdata=[mode1n;mode3n;mode3f];


[nr,nc]=size(mode1n);
[nr1,nc1]=size(Tdata);
[~,mu,sigma]=zscore(mode1n);
mu=repmat(mu,nr1,1);
sigma=repmat(sigma,nr1,1);

% [nr,nc]=size(mode3n);
% [nr1,nc1]=size(Tdata);
% [~,mu,sigma]=zscore(mode3n);
% mu=repmat(mu,nr1,1);
% sigma=repmat(sigma,nr1,1);

dataP1=(Tdata-mu)./sigma;
[T1,P1] = pca(dataP1,10);

% label=zeros(3000,1);
% label(1:2000,:)=classlabel;
% label(2001:3000,:)=3;


label=zeros(3000,1);
% label(1:2000,:)=classlabel1;
% label(2001:4000,:)=classlabel2;
label(1:1000,:)=1;
label(1001:3000,:)=classlabel+1;

[Trans,Z]=LFDA(T1',label,3);
scatter3(Z(1,1:1000),Z(2,1:1000),Z(3,1:1000),'og');
hold on;
% scatter3(Z(1,1001:2000),Z(2,1001:2000),Z(3,1001:2000),'xr');
% scatter3(Z(1,2001:3000),Z(2,2001:3000),Z(3,2001:3000),'^b');
% scatter3(Z(1,3001:4000),Z(2,3001:4000),Z(3,3001:4000),'*m');
scatter3(Z(1,1001:2000),Z(2,1001:2000),Z(3,1001:2000),'^b');
scatter3(Z(1,2001:3000),Z(2,2001:3000),Z(3,2001:3000),'*m');
axis tight 
xlabel('1st Fisher');
ylabel('2nd Fisher');
zlabel('3rd Fisher');
% legend('Base Mode Normal','Base Mode Faulty','New Mode Normal','New Mode Faulty');
legend('Base Mode Normal','New Mode Normal','New Mode Faulty');
% axis([-100 1200 -200 300 -400 0]);
hold off;
