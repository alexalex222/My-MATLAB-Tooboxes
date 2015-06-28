% The split merge inference routine
function [free_energy, posterior] = split_merge(data,posterior,prior,options)
    c_max = 5;
    % first update the current posterior with the new data
    [free_energy, posterior, q_z] = iterate_posterior(data,posterior,prior,options);
    fe_improved = true;
    K = size(q_z.singlets,2) - 1;
    
    iter = 1;
    while fe_improved
        disp('Choosing split or merge...');
        [split_delta_free_energy,member_singlets,member_clumps,sub_partition] = ...
                choose_split_merge(data,posterior,prior,options,q_z);
        
        [foo,split_ind] = sort(split_delta_free_energy,'descend');
              
        for ii = 1:min(c_max,length(split_ind)) 
            split_c = split_ind(ii);
            
            if isempty(sub_partition{split_c})
                fe_improved = false;
                continue;
            end
            
            new_q_z = q_z;
            new_q_z.singlets(member_singlets{split_c},:) = 0;
            new_q_z.singlets(sub_partition{split_c}.member_s1,split_c) = 1;
            new_q_z.singlets(sub_partition{split_c}.member_s2,K+1) = 1;
            new_q_z.singlets(:,K+2) = 0;
        
            new_q_z.clumps(member_clumps{split_c},:) = 0;
        	new_q_z.clumps(sub_partition{split_c}.member_c1,split_c) = 1;
        	new_q_z.clumps(sub_partition{split_c}.member_c2,end) = 1;
        	new_q_z.clumps(:,K+2) = 0;
            
            new_posterior = eval_posterior(data, new_q_z, prior, options); 
            
            [new_free_energy, new_posterior, new_q_z] = iterate_posterior(data,new_posterior,prior,options);
           
            if free_energy_improved(free_energy, new_free_energy, options) 
                free_energy = new_free_energy;
                posterior = new_posterior;
                fe_improved = true;
                K = size(new_q_z.singlets,2) - 1;
                q_z = new_q_z;
                break;
            elseif ii == min(c_max,length(split_ind))
                fe_improved = false;
            end
            
        end
        
    end
    
end
%%
function [split_delta_free_energy,member_singlets,member_clumps,sub_partition] = ...
                            choose_split_merge(data, posterior, prior, options, old_q_z)
    
    K = size(posterior.m,2) - 1;
    c_max = K;
    %merge_candidates = merge_heuristic(data,old_q_z,K);
    %split_candidates = split_heuristic(data,old_q_z,K,posterior,options);
    split_candidates = 1:K;
    
    % book keeping variables
    sub_partition = cell(1,K);
    member_singlets = cell(1,K);
    member_clumps   = cell(1,K);
    
    N_k = zeros(1,K);
    
    % begin by hard assigning to partitions, this is only done to quickly evaluate
    % potential split and merge proposals using local change in free energy; it doesn't change the model
    [foo,top_clusters] = max(old_q_z.singlets,[],2);
    if find(data.Nc) %if there are clumps (i.e. second and later epochs)
        [foo,top_clusters_cl] = max(old_q_z.clumps,[],2);
    end
    for c = 1:K
        member_singlets{c} = find(top_clusters==c);
        N_k(c) = length(member_singlets{c});
        if find(data.Nc) %if there are clumps
            member_clumps{c} = find(top_clusters_cl==c);
            N_k(c) = N_k(c) + sum(data.Nc(member_clumps{c}));
        end
    end
    
    split_delta_free_energy = -inf*ones(1,c_max);
    %merge_delta_free_energy = inf*ones(1,c_max);
    
    % Rank split proposals
    for ii = 1:min(c_max,length(split_candidates))
        
        c = split_candidates(ii);
        
        disp(['Splitting ' num2str(c)]);
        if isempty(member_singlets{c}) && isempty(member_clumps{c})
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
            continue
        end
        display('Proposing component split.');
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
            continue;
        end

        %compute the change in free energy
        split_delta_free_energy(ii) = eval_delta_fe(subdata,sub_posterior,pre_posterior,...
            N_k,new_N_k,options,prior,c);
        %insert_indices = [c (K+1) (K+2)];
        %split_delta_free_energy(ii) = delta_partial_fe(subdata, sub_posterior, pre_posterior, ...
                                       % prior, N_k, new_N_k, insert_indices, options);
                                       
    end
    % Not using merges right now, they don't seem to help much.  
    
    % Rank merge proposals
    %{
    if K > 1
        for tt = 1:min(c_max,size(merge_candidates,2))
            c_i = merge_candidates(1,tt);  c_j = merge_candidates(2,tt);
            disp(['Merging ' num2str(c_i) ' and ' num2str(c_j) '...']);

            merge_singlets = [member_singlets{c_i}; member_singlets{c_j}];
            merge_clumps   = [member_clumps{c_i}; member_clumps{c_j}];
            subdata.singlets = data.singlets(:,merge_singlets);
            subdata.sum_x    = data.sum_x(:,merge_clumps);
            if options.DIAG
                subdata.sum_xx = data.sum_xx(:,merge_clumps);
            else
                subdata.sum_xx = data.sum_xx(:,:,merge_clumps);
            end
            subdata.Nc = data.Nc(merge_clumps);

            % the merged component
            merge_q_z.singlets = zeros(size(merge_singlets,1),2);
            merge_q_z.singlets(:,1) = 1;
            merge_q_z.clumps   = zeros(size(merge_clumps,1),2);
            merge_q_z.clumps(:,1) = 1;

            merge_posterior = eval_posterior(subdata,merge_q_z,prior,options);

            % the split partitions 
            pre_q_z.singlets = zeros(size(merge_singlets,1),3);
            pre_q_z.singlets(1:size(member_singlets{c_i},1),1) = 1;
            pre_q_z.singlets(size(member_singlets{c_i},1)+1:end,2) = 1;
            pre_q_z.clumps = zeros(size(merge_clumps,1),3);
            pre_q_z.clumps(1:size(member_clumps{c_i},1),1) = 1;
            pre_q_z.clumps(size(member_clumps{c_i},1)+1:end,2) = 1;

            pre_posterior = eval_posterior(subdata,pre_q_z,prior,options);

            new_N_k = N_k;
            new_N_k(c_i) = length(merge_singlets) + sum(subdata.Nc);
            new_N_k(c_j) = [];

            merge_delta_free_energy(tt) = eval_delta_fe(subdata, merge_posterior, pre_posterior, ...
                             N_k, new_N_k, options, prior, [c_i c_j]);

        end
    end
    %}
    %{
    [split_free_energy, split_ind] = max(split_delta_free_energy);
    split_c = split_candidates(split_ind);
    %[merge_free_energy, merge_t] = max(merge_delta_free_energy);
    
    %if there's nothing else to do
    %if ~any(~isinf(split_delta_free_energy))&&~any(~isinf(merge_delta_free_energy))
    %    component_ind = [];
    %    return;
    %end
    if ~any(~isinf(split_delta_free_energy))
        component_ind = [];
        return;
    end
    
    %if split_free_energy > merge_free_energy
        component_ind = split_c;
            
        old_q_z.singlets(member_singlets{split_c},:) = 0;
        old_q_z.singlets(sub_partition{split_c}.member_s1,split_c) = 1;
        old_q_z.singlets(sub_partition{split_c}.member_s2,K+1) = 1;
        old_q_z.singlets(:,K+2) = 0;
        
        old_q_z.clumps(member_clumps{split_c},:) = 0;
        old_q_z.clumps(sub_partition{split_c}.member_c1,split_c) = 1;
        old_q_z.clumps(sub_partition{split_c}.member_c2,K+1) = 1;
        old_q_z.clumps(:,K+2) = 0; 
        
    %else
    %    insert_index = merge_candidates(1,merge_t);
    %    delete_index = merge_candidates(2,merge_t);
    %    old_q_z.singlets(member_singlets{insert_index},:) = 0;
    %    old_q_z.singlets(member_singlets{delete_index},:) = 0;
    %    old_q_z.singlets(member_singlets{insert_index},insert_index) = 1;
    %    old_q_z.singlets(member_singlets{delete_index},insert_index) = 1;
    %    old_q_z.singlets(:,delete_index) = [];
        
    %    old_q_z.clumps(member_clumps{insert_index},:) = 0;
    %    old_q_z.clumps(member_clumps{delete_index},:) = 0; 
    %    old_q_z.clumps(member_clumps{insert_index},insert_index) = 1;
    %    old_q_z.clumps(member_clumps{delete_index},insert_index) = 1; 
    %    old_q_z.clumps(:,delete_index) = [];
        
    %    component_ind = [insert_index delete_index];
    %end
%    delta_free_energy = max(split_free_energy,merge_free_energy);
    delta_free_energy = split_free_energy;
    posterior = eval_posterior(data, old_q_z, prior, options);
    %}
end
%% Split heuristic (Ueda et al)
function split_candidates = split_heuristic(data,q_z,K,posterior,options)
    N_k = sum(q_z.singlets,1) + data.Nc*q_z.clumps;
    split_score = zeros(1,K);
    
    D = options.D;
    clump_ind = (data.Nc>0);
    num_clumps = length(find(data.Nc));
    num_singlet = size(data.singlets,2);

    if num_clumps > 0
        sum_x = repmat(reshape(data.sum_x(:,clump_ind),D,1,num_clumps),[1,D,1]);
        Na    = repmat(reshape(data.Nc(clump_ind),1,1,num_clumps),[D,D,1]);
    end
    psi_sum = sum(psi( repmat(posterior.eta+1,D,1) - repmat([1:D]',1,K+1)*0.5 ), 1);
    
    for k = 1:K
        nz_singlets = q_z.singlets(:,k)>0;
        split_score(k) = sum( (q_z.singlets(nz_singlets,k)/N_k(k)) ...
                              .*log(q_z.singlets(nz_singlets,k)/N_k(k)) );
        
        nz_clumps = (data.Nc'.*q_z.clumps(:,k))>0;
        split_score(k) = split_score(k) + sum(data.Nc(nz_clumps)'.*q_z.clumps(nz_clumps,k)/N_k(k) ...
                                                                   .*log(q_z.clumps(nz_clumps,k)/N_k(k)));
                                                               
        E_log_p_of_x_singlets = zeros(1,num_singlet);
        E_log_p_of_x_clumps   = zeros(1,num_clumps);
        if options.DIAG
            precision  = 0.5*diag(posterior.inv_B(:,k))*posterior.eta(k);
            log_par    = -0.5*D*log(pi) - 0.5*sum(log(posterior.B(:,k))) ...
                         + 0.5*psi_sum(k) - 0.5*D/(posterior.xi(k));
            helper_ELPX(E_log_p_of_x_singlets, diag(precision), data.singlets, posterior.m(:,k), log_par);
            if num_clumps > 0
                helper_ELPX_clump(E_log_p_of_x_clumps,data.sum_xx(:,clump_ind),...
                              data.sum_x(:,clump_ind),diag(precision),posterior.m(:,k), log_par, data.Nc(clump_ind));
            end
        else
            precision = 0.5*posterior.inv_B(:,:,k)*posterior.eta(k);
            log_par = -0.5*D*log(pi)  - 0.5*detln(posterior.B(:,:,k)) ...
                      + 0.5*psi_sum(k) - 0.5*D/(posterior.xi(k));
             
            d = data.singlets - repmat(posterior.m(:,k),1,num_singlet);
            E_log_p_of_x_singlets = -sum(d.*(precision*d),1) + log_par;
            
            %helper_ELPX_full(E_log_p_of_x_singlets, precision, data.singlets, posterior.m(:,k), log_par);
            
            if num_clumps > 0
                t2 = sum_x.*repmat(posterior.m(:,k)',[D,1,num_clumps]);
                term_dependent_on_n = (data.sum_xx(:,:,clump_ind) - t2 - permute(t2,[2,1,3]))./Na + ...
                    repmat(posterior.m(:,k)*posterior.m(:,k)',[1,1,num_clumps]);

                E_log_p_of_x_clumps = ...
                    -squeeze(sum(sum(repmat(precision,[1,1,num_clumps]).*term_dependent_on_n,2),1)) + log_par;
            end
        end
                                                               
        split_score(k) = split_score(k) ...
                       - sum((q_z.singlets(:,k)/N_k(k)).*E_log_p_of_x_singlets');
        if num_clumps > 0
            split_score(k) = split_score(k) ...
                       - sum(data.Nc(clump_ind)'.*(q_z.clumps(clump_ind,k)/N_k(k)).*E_log_p_of_x_clumps);
        end
                                                                               
    end
    
    [foo,split_candidates] = sort(split_score,'descend');
    
end
%% Sort all possible candidate pairs according to the merge criterion in Ueda et al
function candidates = merge_heuristic(data,q_z,K)
    merge_crit = q_z.singlets(:,1:K)'*q_z.singlets(:,1:K) + (data.Nc(ones(1,K),:)'.*q_z.clumps(:,1:K))'*q_z.clumps(:,1:K);
    q_norm = sqrt(sum(q_z.singlets(:,1:K).^2,1) + sum(data.Nc(ones(1,K),:)'.*(q_z.clumps(:,1:K).^2),1));
    merge_crit = merge_crit./(q_norm'*q_norm);
    
    if ~isempty(find( (merge_crit).*tril(ones(K),-1) ))
        [all_pairs(1,:),all_pairs(2,:),values] = find( (merge_crit).*tril(ones(K),-1) );
        [foo,ind] = sort(values,'descend');
        candidates = all_pairs(:,ind);
    else
        candidates = [];
    end 
end