function posterior = eval_posterior(data,q_z,prior,options)
    % data magnification factor, used during compression phase only
    if isfield(options,'mag_factor')
        mag_factor = options.mag_factor;
    else
        mag_factor = 1;
    end

    threshold_for_N = 1.0e-200;
    K = size(q_z.singlets,2);
    D = options.D;
    n_items = size(data.singlets,2) + size(data.sum_x,2);
    
    true_N_k = sum(q_z.singlets,1);
    if (~isempty(data.Nc))&find(data.Nc>0)
        true_N_k = true_N_k + data.Nc*q_z.clumps;
    end
    
    q_z.singlets(:,end) = 0;
    q_z.clumps(:,end) = 0;
    
    N_k = mag_factor*true_N_k;

    sum_x = mag_factor*(data.singlets*q_z.singlets + data.sum_x*q_z.clumps);
    
    I = find(N_k>threshold_for_N);
    inv_N_k = zeros(1,K);
    inv_N_k(I) = 1./N_k(I);
    posterior.eta = prior.eta0 + N_k;
    posterior.xi  = prior.xi0  + N_k;
    
    means = sum_x.*repmat(inv_N_k, D, 1);
    
    if options.DIAG
        posterior.inv_B = zeros(D,K);
        posterior.B = zeros(D,K);
    else
        posterior.inv_B = zeros(D,D,K);
        posterior.B = zeros(D,D,K);
    end
        
    for c = 1:K
        v0 = means(:,c) - prior.m0;
        
        if options.DIAG
            if ~isempty(q_z.clumps)
                S = mag_factor*data.sum_xx*q_z.clumps(:,c) ...
                    + mag_factor*helper_B_diag(data.singlets,q_z.singlets(:,c)) ...
                    - N_k(c)*means(:,c).^2;
            else
                S = mag_factor*helper_B_diag(data.singlets,q_z.singlets(:,c)) ...
                    - N_k(c)*means(:,c).^2;
            end
            
             posterior.B(:,c) = prior.B0 + S ...
                            + diag(N_k(c)*prior.xi0*v0*v0')/(posterior.xi(c));
            
        else
            if ~isempty(q_z.clumps)
                cq_of_z_c = reshape(q_z.clumps(:,c), 1, 1, size(q_z.clumps,1));
                S = mag_factor*sum(repmat(cq_of_z_c,[D,D,1]).*data.sum_xx,3) ...
                    + mag_factor*helper_B_full(data.singlets,q_z.singlets(:,c)) ...
                    - N_k(c)*means(:,c)*means(:,c)';
            else
                S = mag_factor*helper_B_full(data.singlets,q_z.singlets(:,c)) - N_k(c)*means(:,c)*means(:,c)';
            end
            posterior.B(:,:,c) = prior.B0 + S ...
                            + N_k(c)*prior.xi0*v0*v0'/(posterior.xi(c)); 
            
        end
        
    end
    
    %also need the inverse of B
    if options.DIAG
        posterior.inv_B = 1./posterior.B;
    else
        for c = 1:K
            posterior.inv_B(:,:,c) = inv(posterior.B(:,:,c));
        end
    end
    
    posterior.m = (sum_x + repmat(prior.xi0*prior.m0,1,K)) ...
                        ./ repmat(N_k+prior.xi0,D,1);
    
    % update the the stick breaking beta distributions
    posterior.gamma = zeros(2,K);
    posterior.gamma(1,:) = 1 + N_k;
    posterior.gamma(2,:) = prior.alpha + sum(N_k) - cumsum(N_k);
    
    posterior.N_k = N_k;
    posterior.true_N_k = true_N_k;
    
end