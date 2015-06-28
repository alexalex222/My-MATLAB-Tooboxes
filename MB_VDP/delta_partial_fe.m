function [delta_free_energy] = delta_partial_fe(data, after_posterior, before_posterior, ...
                                                      hp_prior, before_Nc, after_Nc, insert_indices, opts)
    
    if isfield('opts','mag_factor')
        mag_factor = opts.mag_factor;
    else
        mag_factor = 1;
    end
    
    before_Nc = [before_Nc 0];
    after_Nc = [after_Nc 0];
    
    after_fc    = eval_component_KL(after_posterior, hp_prior, opts); 
    before_fc   = eval_component_KL(before_posterior, hp_prior, opts);
    %after_log_lambda = mk_log_lambda(data, after_posterior, hp_prior,opts); 
    
    before_Nc   = mag_factor*before_Nc;               
    
    before_gamma =  zeros(2,size(before_Nc,2));
    before_gamma(1,:) = 1 + before_Nc;
    before_gamma(2,:) = hp_prior.alpha + sum(before_Nc) - cumsum(before_Nc,2);
    
    before_posterior.gamma(:,1) = before_gamma(:,insert_indices(1));
    before_log_lambda = mk_log_lambda(data, before_posterior, hp_prior, opts);
    
    
    before_E_log_p_of_V = gammaln(sum(before_gamma, 1)) ...
                 - gammaln(1+hp_prior.alpha) ...
                 - sum(gammaln(before_gamma), 1) ...
                 + gammaln(hp_prior.alpha) ...
                 + ((before_gamma(1,:)-1) ...
                    .*(psi(before_gamma(1,:))-psi(sum(before_gamma,1)))) ...
                 + ((before_gamma(2,:)-hp_prior.alpha) ...
                    .*(psi(before_gamma(2,:))-psi(sum(before_gamma,1))));

    after_Nc = mag_factor*after_Nc;
                
    after_gamma = zeros(2,size(after_Nc,2));
    after_gamma(1,:) = 1 + after_Nc;
    after_gamma(2,:) = hp_prior.alpha + sum(after_Nc) - cumsum(after_Nc,2);
    
    %NO RENORMALIZE!
    %after_posterior.gamma = after_gamma(:,insert_indices);    
    
    %RENORMALIZING NOW
    after_posterior.gamma = after_gamma(:,insert_indices);    
    
    after_log_lambda = mk_log_lambda(data, after_posterior, hp_prior,opts);
    
    after_E_log_p_of_V = gammaln(sum(after_gamma, 1)) ...
                 - gammaln(1+hp_prior.alpha) ...
                 - sum(gammaln(after_gamma), 1) ...
                 + gammaln(hp_prior.alpha) ...
                 + ((after_gamma(1,:)-1) ...
                    .*(psi(after_gamma(1,:))-psi(sum(after_gamma,1)))) ...
                 + ((after_gamma(2,:)-hp_prior.alpha) ...
                    .*(psi(after_gamma(2,:))-psi(sum(after_gamma,1))));
            
    
                
    %extra_term = sum(E_log_p_of_V);
    Ns  = size(data.singlets,2);
    
    %log_sum_lambda = log_sum_exp(log_lambda,2);
    before_log_sum_lambda = log_sum_exp(before_log_lambda,2);
    after_log_sum_lambda  = log_sum_exp(after_log_lambda,2);
    
    if ~isempty(data.Nc)
    %    free_energy = extra_term + sum(fc) - sum(log_sum_lambda(1:Ns)) - data.N*log_sum_lambda(Ns+1:end);
        delta_free_energy  = sum(after_E_log_p_of_V) - sum(before_E_log_p_of_V) ...
                       + sum(after_fc) - sum(before_fc) ...
                       + mag_factor*(-sum(after_log_sum_lambda(1:Ns)) + sum(before_log_sum_lambda(1:Ns)) ...
                                     - data.Nc*after_log_sum_lambda(Ns+1:end) + data.Nc*before_log_sum_lambda(Ns+1:end));
    else
    %    free_energy = extra_term + sum(fc) - sum(log_sum_lambda(1:Ns));
        delta_free_energy  = sum(after_E_log_p_of_V) - sum(before_E_log_p_of_V) + ...
                       + sum(after_fc) - sum(before_fc) ...
                       + mag_factor*(-sum(after_log_sum_lambda(1:Ns)) + sum(before_log_sum_lambda(1:Ns)));
    end
end

function log_lambda = mk_log_lambda(data, hp_posterior, hp_prior,opts)

    if abs(hp_posterior.gamma(2,end) - hp_prior.alpha) > 1.0e-5
        hp_posterior.gamma(2,end)
        hp_prior.alpha
        diff = hp_prior.alpha - hp_posterior.gamma(2,end)
        error('must be alpha')
    end
    
    D = size(data.singlets,1);
    N = size(data.singlets,2) + size(data.sum_x,2);
    
    K = size(hp_posterior.eta, 2);

    psi_sum = sum(psi( repmat(hp_posterior.eta+1,D,1) - repmat([1:D]',1,K)*0.5 ), 1); 
    log_lambda = zeros(N,K);
    
    %sum_x = repmat(reshape(data.sum_x,D,1,N),[1,D,1]);
    %Na = repmat(reshape(data.N,1,1,N),[D,D,1]);
    
    clump_ind   = (data.Nc>0);
    %singlet_ind = (data.N == 1);
    %num_singlet = length(find(singlet_ind));
    %num_clumps  = length(find(clump_ind));
    num_singlet = size(data.singlets,2);
    num_clumps   = length(find(data.Nc > 0));
    %num_clumps   = size(data.sum_x,2);
    %xx_ind      = find(clump_ind) - max(find(singlet_ind));
    
    if num_clumps > 0
        sum_x = repmat(reshape(data.sum_x(:,clump_ind),D,1,num_clumps),[1,D,1]);
        Na    = repmat(reshape(data.Nc(clump_ind),1,1,num_clumps),[D,D,1]);
    end
    
    %sum_x = repmat(reshape(data.sum_x(:,clump_ind),D,1,num_clumps),[1,D,1]);
    %Na    = repmat(reshape(data.N(clump_ind),1,1,num_clumps),[D,D,1]);
    
    
    %log_lambda = zeros(N,K);
    for c=1:K  
        E_log_p_of_z_given_other_z_c = psi(hp_posterior.gamma(1,c)) ...
                                     - psi(sum(hp_posterior.gamma(:,c),1)) ...
                                     + sum(psi(hp_posterior.gamma(2,[1:c-1])) ...
                                           - psi(sum(hp_posterior.gamma(:,[1:c-1]),1)), 2);

        Precision = 0.5*hp_posterior.inv_B(:,:,c)*hp_posterior.eta(c);
        pre_E_log_p_of_x = - 0.5*D*log(pi) - 0.5*detln(hp_posterior.B(:,:,c)) ...
                            + 0.5*psi_sum(c) - 0.5*D/(hp_posterior.xi(c));
        
        %E_log_p_of_x  = zeros(1,N);               
        if num_clumps > 0

            t2 = sum_x.*repmat(hp_posterior.m(:,c)',[D,1,num_clumps]);
        
            term_dependent_on_n = (data.sum_xx(:,:,clump_ind) - t2 - permute(t2,[2,1,3]))./Na + ...
                                  repmat(hp_posterior.m(:,c)*hp_posterior.m(:,c)',[1,1,num_clumps]);
        
            d = data.singlets - repmat(hp_posterior.m(:,c),1,num_singlet);
            clump_part   =  -sum(sum(repmat(Precision,[1,1,num_clumps]).*term_dependent_on_n,2),1);
            singlet_part = -sum(d.*(Precision*d),1);
            %E_log_p_of_x  = zeros(1,N);
            E_log_p_of_x.singlets = singlet_part + pre_E_log_p_of_x;
            E_log_p_of_x.clumps  = clump_part + pre_E_log_p_of_x;
            log_lambda(num_singlet+find(clump_ind),c) = E_log_p_of_x.clumps + E_log_p_of_z_given_other_z_c;
        else
            %d = data.sum_x - repmat(hp_posterior.m(:,c),1,N);
            d = data.singlets - repmat(hp_posterior.m(:,c),1,num_singlet);
            E_log_p_of_x.singlets = -sum(d.*(Precision*d),1) + pre_E_log_p_of_x;
             
        end
        
        log_lambda(1:num_singlet,c) = E_log_p_of_x.singlets + E_log_p_of_z_given_other_z_c;
        
    end
    
    log_lambda(:,end) = log_lambda(:,end) - log(1- exp(psi(hp_prior.alpha) - psi(1+hp_prior.alpha)));
    
end


