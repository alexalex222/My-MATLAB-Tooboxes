function component_KL = eval_component_KL(posterior,prior,options)
    D = size(posterior.m, 1);
    K = size(posterior.eta, 2);
    log_det_B = zeros(1,K);
    term_eta = zeros(2,K);
    
    for c=1:K
        d = posterior.m(:,c) - prior.m0;
        if options.DIAG
            log_det_B(c)  = sum(log(posterior.B(:,c)));
            term_eta(1,c) = sum(posterior.inv_B(:,c).*diag(prior.xi0*d*d'),1);
            term_eta(2,c) = sum(posterior.inv_B(:,c).*prior.B0) - D;
        else
            log_det_B(c) = detln(posterior.B(:,:,c));
            term_eta(1,c) = sum(sum(posterior.inv_B(:,:,c).*(prior.xi0*d*d'),1),2);
            term_eta(2,c) = sum(sum(posterior.inv_B(:,:,c).*prior.B0,1),2) - D;
        end
    end
    E_log_q_p_mean = 0.5*D*(prior.xi0./posterior.xi ...
                            - log(prior.xi0./posterior.xi) - 1) ...
                   + 0.5*(posterior.eta).* term_eta(1,:);

    psi_sum = sum(psi( (repmat(posterior.eta+1,D,1) - repmat([1:D]',1,K))*0.5 ), 1); 
    
    if options.DIAG
        log_det_B0 = sum(log(prior.B0));
    else
        log_det_B0 = detln(prior.B0);
    end
    
    E_log_q_p_cov = 0.5*prior.eta0*(log_det_B-log_det_B0) ...
                  + 0.5*posterior.N_k.*psi_sum + 0.5*(posterior.eta).* term_eta(2,:) ...
                  + gamma_multivariate_ln(prior.eta0*0.5,D) ...
                  - gamma_multivariate_ln(posterior.eta*0.5,D);
              
   component_KL = E_log_q_p_mean + E_log_q_p_cov; 
end