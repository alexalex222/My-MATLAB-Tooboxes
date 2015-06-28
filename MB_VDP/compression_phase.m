function data = compression_phase(data,prior,posterior,model_q_z,options,T)
    K = size(posterior.m,2) - 1;
    
    display('Compression phase.');
    options.mag_factor = options.N/T;
    
    % book keeping variables
    sub_partition = cell(1,K);
    member_singlets = cell(1,K);
    member_clumps   = cell(1,K);
    evaluated = zeros(1,K);
    delta_free_energy = -inf*ones(1,K);
    N_k = zeros(1,K); 
    
    % begin by hard assigning to partitions
    [foo,top_clusters] = max(model_q_z.singlets,[],2);
    if find(data.Nc) %if there are clumps (i.e. second and later epochs)
        [foo,top_clusters_cl] = max(model_q_z.clumps,[],2);
    end
    for c = 1:K
        member_singlets{c} = find(top_clusters==c);
        N_k(c) = length(member_singlets{c});
        if find(data.Nc) %if there are clumps
            member_clumps{c} = find(top_clusters_cl==c);
            N_k(c) = N_k(c) + sum(data.Nc(member_clumps{c}));
        end
    end
    
    THRESHOLD = (options.D+3)/2; % don't bother to split clumps below this size (points will be singlets)
        
    while memory_cost(N_k,options,THRESHOLD) < options.M
        for c = 1:K
            if N_k(c) < THRESHOLD % no need to split further
                evaluated(c) = true;
            end
            
            if ~evaluated(c)
                if isempty(member_singlets{c}) && isempty(member_clumps{c})
                    evaluated(c) = true;
                    continue;
                end
                
                subdata.singlets = data.singlets(:,member_singlets{c});
                subdata.sum_x    = data.sum_x(:,member_clumps{c});
                if options.DIAG
                    subdata.sum_xx   = data.sum_xx(:,member_clumps{c});
                else
                    subdata.sum_xx   = data.sum_xx(:,:,member_clumps{c});
                end
                subdata.Nc       = data.Nc(member_clumps{c});
                
                %construct partition parameters (before split)
                pre_q_z.singlets = zeros(size(member_singlets{c},1),2);
                pre_q_z.singlets(:,1) = 1;
                pre_q_z.clumps = zeros(size(member_clumps{c},1),2);
                pre_q_z.clumps(:,1) = 1;
                pre_posterior = eval_posterior(subdata,pre_q_z,prior,options);
                
                % propose split along principal componenent
                sub_q_z = split(subdata,1,pre_q_z,pre_posterior,options);
                % if this fails to produce a split
                if ~isempty(find((sum(sub_q_z.singlets(:,1:end-1),1)+sum(sub_q_z.clumps((subdata.Nc>0),1:end-1),1))<1.0e-10,1))
                    evaluated(c) = true;
                    continue
                end
                display('Proposing partition split.');
                % refine the partition with local updates
                sub_posterior = eval_posterior(subdata,sub_q_z,prior,options);
                [sub_fe, sub_posterior, sub_q_z] = iterate_posterior(subdata,sub_posterior,prior,options,10,false,true);
                
                %insert_indices = [c (K+1) (K+2)];
                new_N_k = N_k;
                                
                %store subpartitions
                [foo,top_sub_par] = max(sub_q_z.singlets,[],2);
                sub_partition{c}.member_s1 = member_singlets{c}(top_sub_par==1);
                new_N_k(c) = length(sub_partition{c}.member_s1);
                sub_partition{c}.member_s2 = member_singlets{c}(top_sub_par==2);
                new_N_k(K+1) = length(sub_partition{c}.member_s2);
                
                if find(subdata.Nc)
                    [foo,top_sub_par_cl] = max(sub_q_z.clumps,[],2);
                    sub_partition{c}.member_c1 = member_clumps{c}(top_sub_par_cl==1);
                    new_N_k(c) = new_N_k(c) + sum(data.Nc(sub_partition{c}.member_c1));
                    sub_partition{c}.member_c2 = member_clumps{c}(top_sub_par_cl==2);
                    new_N_k(K+1) = new_N_k(K+1) + sum(data.Nc(sub_partition{c}.member_c2));
                else
                    sub_partition{c}.member_c1 = [];
                    sub_partition{c}.member_c2 = [];
                end
                % if the partition didn't split
                if ((isempty(sub_partition{c}.member_s1) && isempty(sub_partition{c}.member_c1)) || ...
                    (isempty(sub_partition{c}.member_s2) && isempty(sub_partition{c}.member_c2)))
                    evaluated(c) = true;
                    continue;
                end
                
                %compute the change in free energy                
                delta_free_energy(c) = eval_delta_fe(subdata,sub_posterior,pre_posterior,...
                                                     N_k,new_N_k,options,prior,c);
                %insert_indices = [c (K+1) (K+2)];
                
                %delta_free_energy(c) = delta_partial_fe(subdata, sub_posterior, pre_posterior, ...
                %                                      prior, N_k, new_N_k, insert_indices, options);
                
                
                evaluated(c) = true;
            end
             
        end
        
        [free_energy,c] = max(delta_free_energy);
        if isinf(delta_free_energy)
            break;
        end
        
        member_singlets{c} = sub_partition{c}.member_s1;
        member_clumps{c}   = sub_partition{c}.member_c1;
        member_singlets{K+1} = sub_partition{c}.member_s2;
        member_clumps{K+1} = sub_partition{c}.member_c2;
        
        N_k(c) = length(member_singlets{c}) + sum(data.Nc(member_clumps{c}));
        N_k(K+1) = length(member_singlets{K+1}) + sum(data.Nc(member_clumps{K+1}));
        
        evaluated(c) = false;
        evaluated(K+1) = false;
        delta_free_energy([c K+1]) = -inf;
        K = K + 1;
    end
    
    % summarize the new partitions with their sufficient statistics
    the_clumps = (N_k >= THRESHOLD);
    the_singlets = (N_k > 0) & (N_k < THRESHOLD);
    
    temp_x = data.sum_x; temp_xx = data.sum_xx; temp_Nc = data.Nc;
    
    for cc = find(the_clumps)
        data.sum_x(:,cc) = sum(data.singlets(:,member_singlets{cc}),2) ...
                         + sum(temp_x(:,member_clumps{cc}),2);
        
        data.Nc(cc) = length(member_singlets{cc}) + sum(temp_Nc(member_clumps{cc}));
        
        if options.DIAG
            data.sum_xx(:,cc) = sum(data.singlets(:,member_singlets{cc}).^2,2) ...
                              + sum(temp_xx(:,member_clumps{cc}),2);
        else
            if ~isempty(member_singlets{cc})
                data.sum_xx(:,:,cc) = ...
                    data.singlets(:,member_singlets{cc})*data.singlets(:,member_singlets{cc})';
            else
                data.sum_xx(:,:,cc) = zeros(options.D,options.D,1);
            end
            if ~isempty(member_clumps{cc})
                data.sum_xx(:,:,cc) = data.sum_xx(:,:,cc) + sum(temp_xx(:,:,member_clumps{cc}),3);
            end 
        end
    end
    
    data.sum_x = data.sum_x(:,the_clumps);
    if options.DIAG
        data.sum_xx = data.sum_xx(:,the_clumps);
    else
        data.sum_xx = data.sum_xx(:,:,the_clumps);
    end
    data.Nc = data.Nc(the_clumps);
    
    data.singlets = data.singlets(:,vertcat(member_singlets{the_singlets}));
    
end

function cost = memory_cost(N_k,options,THRESHOLD)
    if options.DIAG
        CLUMP_COST = 2;
    else
        CLUMP_COST = (options.D+3)/2;
    end
    
    num_clumps = length(find(N_k>=THRESHOLD));
    num_singlets = sum(N_k(N_k<THRESHOLD));
    
    cost = ceil(num_clumps*CLUMP_COST + num_singlets);
end
