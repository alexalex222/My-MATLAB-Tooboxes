function S_sk = eval_S_sk(data,posterior,prior,options,relevant_components)
    %K = size(posterior.m,2);
    
    if nargin < 5
        K = size(posterior.m,2);
    	relevant_components = 1:K;
    else
        K = length(relevant_components);
    end
    D = options.D;
    n_items = size(data.singlets,2) + size(data.sum_x,2); % # of singlets + # of clumps
    
    S_sk = zeros(n_items, K);
    
    clump_ind = (data.Nc>0);
    num_clumps = length(find(data.Nc));
    num_singlet = size(data.singlets,2);
    
    if num_clumps > 0
        sum_x = repmat(reshape(data.sum_x(:,clump_ind),D,1,num_clumps),[1,D,1]);
        Na    = repmat(reshape(data.Nc(clump_ind),1,1,num_clumps),[D,D,1]);
    end
    
    %psi_sum = sum(psi( repmat(posterior.eta+1,D,1) - repmat([1:D]',1,K)*0.5 ), 1);
    %psi_sum = sum(psi( repmat(posterior.eta+1,D,1) - repmat([1:D]',1,length(posterior.eta))*0.5 ), 1);
    psi_sum = sum(psi( (repmat(posterior.eta+1,D,1) - repmat([1:D]',1,length(posterior.eta)))*0.5 ), 1);
   
    
    for c = 1:K
        g_c = relevant_components(c);
        E_log_p_of_z_given_V = psi(posterior.gamma(1,g_c)) - psi(sum(posterior.gamma(:,g_c),1)) ...
                                     + sum(psi(posterior.gamma(2,[1:g_c-1])) ...
                                           - psi(sum(posterior.gamma(:,[1:g_c-1]),1)), 2);

        E_log_p_of_x_singlets = zeros(1,num_singlet);
       
        if options.DIAG 
            E_log_p_of_x_clumps   = zeros(1,num_clumps);
            precision  = 0.5*diag(posterior.inv_B(:,c))*posterior.eta(c);
            log_par    = -0.5*D*log(pi) - 0.5*sum(log(posterior.B(:,c))) ...
                         + 0.5*psi_sum(c) - 0.5*D/(posterior.xi(c));
            helper_ELPX(E_log_p_of_x_singlets, diag(precision), data.singlets, posterior.m(:,c), log_par);
            
            if num_clumps > 0
                helper_ELPX_clump(E_log_p_of_x_clumps,data.sum_xx(:,clump_ind),...
                              data.sum_x(:,clump_ind),diag(precision),posterior.m(:,c), log_par, data.Nc(clump_ind));
                
                S_sk(num_singlet+find(clump_ind),c) = E_log_p_of_x_clumps + E_log_p_of_z_given_V;
           
            end
        else
            precision = 0.5*posterior.inv_B(:,:,c)*posterior.eta(c);
            log_par = -0.5*D*log(pi)  - 0.5*detln(posterior.B(:,:,c)) ...
                      + 0.5*psi_sum(c) - 0.5*D/(posterior.xi(c));
                  
            %helper_ELPX_full(E_log_p_of_x_singlets, precision, data.singlets, posterior.m(:,c), log_par);
            d = data.singlets - repmat(posterior.m(:,c),1,num_singlet);
            E_log_p_of_x_singlets = -sum(d.*(precision*d),1) + log_par;
            
            if num_clumps > 0
                t2 = sum_x.*repmat(posterior.m(:,c)',[D,1,num_clumps]);
                term_dependent_on_n = (data.sum_xx(:,:,clump_ind) - t2 - permute(t2,[2,1,3]))./Na + ...
                    repmat(posterior.m(:,c)*posterior.m(:,c)',[1,1,num_clumps]);

                S_sk(num_singlet+find(clump_ind),c) = ...
                    -squeeze(sum(sum(repmat(precision,[1,1,num_clumps]).*term_dependent_on_n,2),1)) + log_par ...
                    + E_log_p_of_z_given_V;
            end
        end
        
        S_sk(1:num_singlet,c) = E_log_p_of_x_singlets + E_log_p_of_z_given_V;
                                       
    end
    
    S_sk(:,end) = S_sk(:,end) - log(1- exp(psi(prior.alpha) - psi(1+prior.alpha)));
end