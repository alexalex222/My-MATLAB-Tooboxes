%%
load pcs.mat
first_pcs = pcs(:,1:50);
%%
res_str = 'results';
K = eval([res_str '.K']);
m = eval([res_str '.posterior.m']);
IMAGE_SIZE = 28;
PCA = 1;
num_cols = 10;
num_rows = ceil((K-1)/10);

%%
figure;

IMAGE_SIZE = 28;
num_cols = 10;
num_rows = ceil((K-1)/10);

panel = zeros(IMAGE_SIZE*num_rows,IMAGE_SIZE*num_cols);

for t = 1:K-1
  cluster_mean = m(:,t);
  
  projected_vector = first_pcs*cluster_mean + mv';
  the_image = reshape(projected_vector,IMAGE_SIZE,IMAGE_SIZE)';
  the_image(the_image < 0) = 0;
  the_image = the_image/norm(the_image(:));
  
  %the_image = reshape(cluster_mean,7,7);
  
  current_row = ceil(t/num_cols);
  current_col = mod(t,num_cols);
  if current_col == 0
      current_col = num_cols;
  end
    
  panel((current_row-1)*IMAGE_SIZE+1:current_row*IMAGE_SIZE, ...
        (current_col-1)*IMAGE_SIZE+1:current_col*IMAGE_SIZE) = the_image;
  %imagesc(panel),colormap(gray);
  %current_row,current_col
  
end

imagesc(panel),colormap(gray);

hold on; axis equal; axis off;
for i = 0:num_cols %vertical lines
    plot([i*IMAGE_SIZE+.5 i*IMAGE_SIZE+.5],[+.5 num_rows*IMAGE_SIZE+.5]);
end
for i = 0:num_rows %horizontal lines
    plot([+.5 num_cols*IMAGE_SIZE+0.5],[i*IMAGE_SIZE+.5 i*IMAGE_SIZE+.5]);
end
%% clumps
res_str = 'results';
data = eval([res_str '.data']);
D=50;
m = data.sum_x./repmat(data.Nc,D,1);
[foo,ind]=sort(data.Nc,'descend');
m=m(:,ind);
K=length(data.Nc)+1;
IMAGE_SIZE = 28;
PCA = 1;
num_cols = 10;
num_rows = ceil((K-1)/10);

%% Compare cluster profiles of incremental vs batch
%load results_vdp.mat
%figure(1);
%plot(results.hp_posterior.Nc(1:end-1),'b');
%load gs_3000_3000_19.mat
%hold on;
%plot(results_GS.hp_posterior.Nc(1:end-1),'r');
%legend('Batch','Incremental');
%xlabel('Cluster');
%ylabel('Number of data points');
%title('Cluster Size Profile');
%% Assign digits to cluster using the responsibilities
%{
%[resp_sorted,resp_i] = sort(results_ML.q_of_z,2,'descend');
%figure(11); 
%subplot(1,3,1); hist(resp_sorted(:,1),100); title('First responsibilities');
%subplot(1,3,2); hist(resp_sorted(:,2),100); title('Second responsibilities');
%subplot(1,3,3); hist(resp_i(:,1),[1:results_MLC.K]); title('Cluster populations');
%}
%%



%% Now show samples from each cluster
%{
%CLUST_PER_FIG = 39;
%IMA_PER_CLUST  = 12;
%IMAGE_SIZE = 28;
%figure_index =0;
%for t=1:results_ML.K-1, %% Loop over all clusters that were found
%    if mod(t-1,CLUST_PER_FIG)==0,
%        figure_index = figure_index+1;
%        figure(100+figure_index);
%        ima_clusters = zeros(CLUST_PER_FIG*IMAGE_SIZE,IMA_PER_CLUST*IMAGE_SIZE);
%    end;
%    i_imgs_incluster = find(resp_i(:,1)==t);
    %resp_imgs_incluster = resp_sorted(i_imgs_incluster,1);
    %[resp_imgs_sorted,resp_imgs_sorted_i] = sort(resp_imgs_incluster,'descend');
%    i_imgs_displayed = i_imgs_incluster(randsample(length(i_imgs_incluster),IMA_PER_CLUST));
%    resp_imgs_displayed = resp_sorted(i_imgs_displayed,1);
%    k_row = t-CLUST_PER_FIG*floor((t-1)/CLUST_PER_FIG);
%    fprintf(1,'t=%d, k_row=%d k_col=',t, k_row);
%    for k_col = 1:IMA_PER_CLUST,
%        ima_clusters((k_row-1)*IMAGE_SIZE+[1:IMAGE_SIZE],(k_col-1)*IMAGE_SIZE+[1:IMAGE_SIZE]) = ...
%            reshape(reducedFeatures(i_imgs_displayed,:)*first_pcs' + mv(ones(IMA_PER_CLUST,1),:),IMAGE_SIZE,IMAGE_SIZE)';
%        fprintf(1,'%d ',k_col);
%    end
%    fprintf(1,'\n');
%    imagesc(ima_clusters); colormap(gray);
%    for k_col = 1:IMA_PER_CLUST,
%        text((k_col-1)*IMAGE_SIZE+1,(k_row-1)*IMAGE_SIZE+IMAGE_SIZE,num2str(round(10*resp_imgs_displayed(k_col))),'Color','w');
%    end
%end
%}

