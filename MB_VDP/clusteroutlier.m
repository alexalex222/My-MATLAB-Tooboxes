clear;
clc;
close all;

load plsfilled.mat


temp=plsxfilled(:,[2:36,38:62]);

[nr,nc]=size(temp);

[dataP1,mu,sigma]=zscore(temp);






[T1,P1] = pca(dataP1,60);
results = MB_VDP(T1',11787,0,generate_prior(T1'));
[~,classlabel]=max(results.q_z.singlets,[],2);




%% Scatter plot
figure; %Plot unlabeled data
scatter3(T1(:,1),T1(:,2),T1(:,3));


figure; % Plot labeled data
scatter3dMixture(T1,classlabel);
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5');


%% cluster again

% temp2=plsxfilled(classlabel==1,[2:36,38:39,41:62]);

temp2=plsxfilled(classlabel==1,[2:36,38:39,41:62]);

[dataP2,mu,sigma]=zscore(temp2);

[T2,P2] = pca(dataP2,60);
results = MB_VDP(T2',9937,0,generate_prior(T2'));
[~,classlabel2]=max(results.q_z.singlets,[],2);


%% plot again

figure; % Plot labeled data
scatter3dMixture(T2,classlabel2);
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');