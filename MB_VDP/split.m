function new_q_z = split(data,c,q_z,posterior,options)
    if c > size(posterior.B,2) - 1
        cluster_mean = posterior.m(:,1);
        if options.DIAG
            cluster_cov  = posterior.B(:,1)/posterior.eta(1);
        else
            cluster_cov  = posterior.B(:,:,1)/posterior.eta(1);
        end
    else
        cluster_mean = posterior.m(:,c);
        if options.DIAG
            cluster_cov  = posterior.B(:,c)/posterior.eta(c);
        else
            cluster_cov  = posterior.B(:,:,c)/posterior.eta(c);
        end
    end
    
    valid_ind = (data.Nc > 0);
    
    if options.DIAG
        dir_singlet = divide_diagonal(data.singlets, cluster_cov, cluster_mean);
        if find(valid_ind)
            dir_clump   = divide_diagonal(data.sum_x(:,valid_ind)./data.Nc(ones(size(data.sum_x,1),1),valid_ind), ...
                                  cluster_cov, cluster_mean);
        else
            dir_clump = [];
        end
    else
        dir_singlet = divide_by_principal_component(data.singlets, cluster_cov, cluster_mean);
        if find(valid_ind)
            dir_clump   = divide_by_principal_component(data.sum_x(:,valid_ind)./data.Nc(ones(size(data.sum_x,1),1),valid_ind), ...
                                  cluster_cov, cluster_mean);
        else
            dir_clump = [];
        end
    end

    q_z_c1_s = zeros(size(q_z.singlets,1),1);
    q_z_c2_s = q_z.singlets(:,c);
    I = find(dir_singlet>=0);
    q_z_c1_s(I) = q_z.singlets(I,c);
    q_z_c2_s(I) = 0;
    
    new_q_z.singlets = zeros(size(q_z.singlets,1), size(q_z.singlets,2)+1);
    new_q_z.singlets(:,[1:end-2 end]) = q_z.singlets;
    new_q_z.singlets(:,c) = q_z_c1_s;
    new_c = size(new_q_z.singlets, 2) - 1;
    new_q_z.singlets(:,new_c) = q_z_c2_s;
    
    q_z_c1_cl = zeros(size(q_z.clumps,1),1);
    q_z_c2_cl = q_z.clumps(:,c);
    I = find(dir_clump>=0);
    q_z_c1_cl(I) = q_z.clumps(I,c);
    q_z_c2_cl(I) = 0;
    
    new_q_z.clumps = zeros(size(q_z.clumps,1), size(q_z.clumps,2)+1);
    new_q_z.clumps(:,[1:end-2 end]) = q_z.clumps;
    new_q_z.clumps(:,c) = q_z_c1_cl;
    new_c = size(new_q_z.clumps, 2) - 1;
    new_q_z.clumps(:,new_c) = q_z_c2_cl;
    
end
%%
function direction = divide_by_principal_component(data, covariance, mean)
    N = size(data, 2);
    [V,D] = eig(covariance);
    [eig_val, principal_component_i] = max(diag(D));
    principal_component = V(:,principal_component_i);
    direction = sum((data - repmat(mean, 1, N)).*repmat(principal_component, 1, N), 1);
end
%%
function direction = divide_diagonal(data,covariance,mean)
    N = size(data, 2);
    [foo,pc_idx] = max(covariance);
    direction = data(pc_idx,:) - mean(pc_idx);
end