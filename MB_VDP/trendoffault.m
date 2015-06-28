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



mode1n=[mode1_normal.x];
mode3n=[mode3_normal.x];
mode1f=[mode1_idv14.x];
mode3f=[mode3_idv13.x];

Tdata=[mode3n;mode3f];

% Tdata=[mode3n;mode3f];


% %% case1 mode1idv1
% figure
% subplot(3,1,1)
% plot(Tdata(:,1));
% ylabel('Trend');
% xlabel('Variable 1');
% subplot(3,1,2);
% plot(Tdata(:,1));
% ylabel('Trend');
% xlabel('Variable 4');
% subplot(3,1,3);
% plot(Tdata(:,18));
% ylabel('Trend');
% xlabel('Variable 18');


% %% case2 mode1idv14
% subplot(2,1,1)
% plot(Tdata(:,9));
% ylabel('Trend');
% xlabel('Variable 9');
% subplot(2,1,2);
% plot(Tdata(:,21));
% ylabel('Trend');
% xlabel('Variable 21');

%% case3 mode3idv13
subplot(3,2,1)
plot(Tdata(:,7));
ylabel('Trend');
xlabel('Variable 7');
subplot(3,2,2)
plot(Tdata(:,11));
ylabel('Trend');
xlabel('Variable 11');
subplot(3,2,3)
plot(Tdata(:,13));
ylabel('Trend');
xlabel('Variable 13');
subplot(3,2,4)
plot(Tdata(:,15));
ylabel('Trend');
xlabel('Variable 15');
subplot(3,2,5)
plot(Tdata(:,16));
ylabel('Trend');
xlabel('Variable 16');
subplot(3,2,6)
plot(Tdata(:,21));
ylabel('Trend');
xlabel('Variable 21');