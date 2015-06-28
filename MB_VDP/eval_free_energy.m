function [free_energy, S_sk] = eval_free_energy(data, posterior, prior, options)
    
    if isfield(options,'mag_factor')
        mag_factor = options.mag_factor;
    else
        mag_factor = 1;
    end
    
    component_KL = eval_component_KL(posterior,prior,options);
    S_sk = eval_S_sk(data,posterior,prior,options);
    
    stick_break_KL = gammaln(sum(posterior.gamma, 1)) ...
                   - gammaln(1+prior.alpha) ...
                   - sum(gammaln(posterior.gamma), 1) ...
                   + gammaln(prior.alpha) ...
                   + ((posterior.gamma(1,:)-1) ...
                     .*(psi(posterior.gamma(1,:))-psi(sum(posterior.gamma,1)))) ...
                   + ((posterior.gamma(2,:)-prior.alpha) ...
                     .*(psi(posterior.gamma(2,:))-psi(sum(posterior.gamma,1))));    

    num_singlets = size(data.singlets,2);
    
    lse_S_sk = helper_logsumexp(S_sk,max(S_sk,[],2));
    
    free_energy = -sum(component_KL) - sum(stick_break_KL) ...
                + mag_factor*sum(lse_S_sk(1:num_singlets));
    
    if ~isempty(data.Nc)
        free_energy = free_energy + mag_factor*data.Nc*lse_S_sk(num_singlets+1:end);
    end

end