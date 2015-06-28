%% Memory Bounded Variational Dirichlet Process 
% X: Each column is a D dimensional data point.  Passing it in as one big matrix
%    will work when the data fits in memory.  For larger problems, modify
%    the main loop to retrieve E data points from disk during each round.
% M: The memory bound, measured as the maximum number of data points that could be
%    stored.  This is used to store clumps that summarize previous learning
%    rounds.  There is an additional buffer to store the current epoch of E
%    data points.
% E: the epoch size (number of points to process at one time).  The first 
%    chunk of data will be M+E in size, and then after that E points will
%    be added at a time.
% prior: prior hyperparameters. The generate_prior.m function can be used to create this.  
% file_name: (optional) file name to write state to after each epoch
function results = MB_VDP(X,M,E,prior,file_name)
    if nargin < 5
        file_name = [];
    end
    
    start_time = clock;
    % relative threshold used to judge free energy convergence
    options.threshold = 1.0e-5;
    % if options.restart=true, then the model is built from scratch (all points
    %    initialized to same cluster) after each learning round.  If options.restart=false
    %    then the previous model estimate is used as a starting point.
    %    This is much faster and typically results in quite similar models.
    options.restart = false; 
    % Plots the data points, model clusters, and clumps (just the first
    % two dimensions), if options.display=true. 
    options.display = false;
    options.M = M;
    options.E = E;
    options.N = size(X,2); %Total size of the data set
    options.D = size(X,1); %Dimensionality of data set
    if size(prior.B0,2) == 1
        options.DIAG = true; %diagonal covariance matrices
    else
        options.DIAG = false; %full covariance matrices
    end
    D = options.D;
    
    %initialize working memory with first data points
    data.singlets = X(:,1:(M+E));
    %clump sufficient statistics
    max_clumps = ceil((2*M)/(D+3));
    data.sum_x = zeros(D,max_clumps);
    if options.DIAG
        data.sum_xx = zeros(D,max_clumps);
    else
        data.sum_xx = zeros(D,D,max_clumps);
    end
    data.Nc = zeros(1,max_clumps);
    
    n_pts_in_epoch = M + E; %the first epoch has more points
    
    % the main loop
    ii = 0;
    T = 0; %the number of datapoints processed so far
    
    while true
        ii = ii + 1;
        display(['Incremental Cycle ' num2str(ii)]);
        
        if options.restart || (ii == 1)
            % initialize all points to same cluster
            q_z.singlets = zeros(size(data.singlets,2),2);
            q_z.singlets(:,1) = 1;
            q_z.clumps   = zeros(size(data.sum_x,2),2);
            q_z.clumps(:,1) = 1;
            % the posterior structure holds model hyperparameters
            posterior = eval_posterior(data,q_z,prior,options);
        end

        [free_energy, posterior] = split_merge(data,posterior,prior,options);
    
        T = T + n_pts_in_epoch;
        
        if T < options.N
            q_z = eval_q_z(data,posterior,prior,options);
            data = compression_phase(data,prior,posterior,q_z,options,T);
            
            % Display the current state on first two dimensions
           if options.display
                currentAxis = [-2 2 -2 2];
                axis(currentAxis);
                axis square;
                cla;
                %plot all data points in blue
                plot(X(1,:),X(2,:),'.b');

                %display discarded points as red
                hold on;
                plot(X(1,1:T),X(2,1:T),'.r');
                hold off;
                %display the model
                displayModel(posterior,'g',options);
                %display the contents of working memory (individual points in
                %  green, compressed statistics as red ellipses
                displayPts(data,'g',options);

                drawnow;
           end
            
            if (options.N - T) >= E
                n_pts_in_epoch = E;
            else
                n_pts_in_epoch = options.N - T;
            end
            
            % get the next epoch of singlets
            data.singlets = [data.singlets X(:,T+1:T+n_pts_in_epoch)];
            
            % write the state to a file
            if (~isempty(file_name))
                results.data = data;
                results.elapsed_time = etime(clock, start_time);
                results.free_energy = free_energy;
                results.prior = prior;
                results.posterior = posterior;
                results.K = size(posterior.m,2);
                results.options = options;
      
                save([file_name '_' num2str(ii) '.mat'],'results');
                
            end    
        else
            break;
        end
    end
    
    results.data = data;
    results.elapsed_time = etime(clock, start_time);
    results.free_energy = free_energy;
    results.prior = prior;
    results.posterior = posterior;
    results.K = length(posterior.eta);
    results.options = options; 
    results.q_z = eval_q_z(data,posterior,prior,options);
    
    % Display 
    if options.display
        currentAxis = [-2 2 -2 2];
        axis(currentAxis);
        axis square;
        cla;
        %plot all data points in blue
        plot(X(1,:),X(2,:),'.b');

        %display discarded points as red
        hold on;
        plot(X(1,1:T),X(2,1:T),'.r');
        hold off;
        %display the model
        displayModel(posterior,'g',options);
        %display the contents of working memory (individual points in
        %  green, compressed statistics as red ellipses
        displayPts(data,'g',options);

        drawnow;
    end
end

%% Display the working memory contents
function displayPts(data,col,options)
    hold on;
    
    plot(data.singlets(1,:),data.singlets(2,:),['.' col]);
    
    for i = 1:size(data.Nc,2)
        if data.Nc(i) > 0
            clump_mean = data.sum_x(:,i)/data.Nc(i);
            if options.DIAG
                clump_cov = diag(data.sum_xx(:,i))/data.Nc(i) ...
                    - diag(diag(clump_mean*clump_mean')) + 0.001*eye(size(data.sum_x,1));
            else
                clump_cov = data.sum_xx(:,:,i)/data.Nc(i) ...
                    - clump_mean*clump_mean' + 0.001*eye(size(data.sum_x,1));
            end
            draw_ellipse(clump_mean(1:2),clump_cov(1:2,1:2),'r');
        end
    end
end
%% Display the current model
function displayModel(posterior,col,options)
    for j = find(posterior.true_N_k >= 1)
        if options.DIAG
            theCov = diag(posterior.B(:,j))/posterior.eta(j);
        else
            theCov = posterior.B(:,:,j)/posterior.eta(j);
        end
        draw_ellipse(posterior.m(1:2,j),theCov(1:2,1:2),col);
    end
end