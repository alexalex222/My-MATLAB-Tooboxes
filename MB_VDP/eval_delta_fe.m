% fast evaluation of the change in free energy when points are hard
% assigned to partitions.  This is used to evaluate split-merge proposals
% and to determine which partitions to split when determining clumps
function delta_fe = eval_delta_fe(subdata,post_posterior,pre_posterior, ...
                                  pre_N_k,post_N_k,options,prior,comps)
    if isfield(options,'mag_factor')
        mag_factor = options.mag_factor;
    else
        mag_factor = 1;
    end
    
    post_component_KL = eval_component_KL(post_posterior, prior, options);
    pre_component_KL = eval_component_KL(pre_posterior, prior, options);
    
    [pre_N_k,pre_idx] = sort(mag_factor*pre_N_k,'descend'); %DP stick weights are size biased
    [post_N_k,post_idx] = sort(mag_factor*post_N_k,'descend');
    
    %foo = zeros(1,length(post_N_k));
    %foo(post_idx) = 1:length(post_N_k);
    % post_idx takes a new component order and returns its old order
    % foo takes an old component order and returns its new order
    
    % reconstruct the stick breaking hyperparameters (pre_gamma could be
    % computed outside of this function to reduce redundancy, but this would
    % make the code more difficult for not much speedup)
    
    pre_N_k = [pre_N_k 0];
    post_N_k = [post_N_k 0];
    
    pre_gamma = zeros(2,length(pre_N_k));
    pre_gamma(1,:) = 1 + pre_N_k;
    pre_gamma(2,:) = prior.alpha + sum(pre_N_k) - cumsum(pre_N_k);
    post_gamma = zeros(2,length(post_N_k));
    post_gamma(1,:) = 1 + post_N_k;
    post_gamma(2,:) = prior.alpha + sum(post_N_k) - cumsum(post_N_k);

    pre_stick_KL = gammaln(sum(pre_gamma, 1)) - gammaln(1+prior.alpha) ...
                        - sum(gammaln(pre_gamma), 1) + gammaln(prior.alpha) ...
                        + ((pre_gamma(1,:)-1).*(psi(pre_gamma(1,:))-psi(sum(pre_gamma,1)))) ...
                        + ((pre_gamma(2,:)-prior.alpha).*(psi(pre_gamma(2,:))-psi(sum(pre_gamma,1))));
    
    post_stick_KL = gammaln(sum(post_gamma, 1)) - gammaln(1+prior.alpha) ...
                        - sum(gammaln(post_gamma), 1) + gammaln(prior.alpha) ...
                        + ((post_gamma(1,:)-1).*(psi(post_gamma(1,:))-psi(sum(post_gamma,1)))) ...
                        + ((post_gamma(2,:)-prior.alpha).*(psi(post_gamma(2,:))-psi(sum(post_gamma,1))));
    
    % put in correct stick breaking weights in order to evaluate S_sk
    pre_posterior.gamma = pre_gamma;
    post_posterior.gamma = post_gamma;
    
    if length(comps) == 1 % split
        %pre_posterior.gamma(:,1) = pre_gamma(:,pre_idx==comps);
        %post_posterior.gamma(:,1) = post_gamma(:,post_idx==comps);
        %post_posterior.gamma(:,2) = post_gamma(:,post_idx==length(post_N_k));
        pre_comps = [find(pre_idx==comps) length(pre_N_k)];
        post_comps = [find(post_idx==comps) find(post_idx==(length(post_N_k)-1)) length(post_N_k)];
        %foo_post_comps = [foo(comps) foo(length(pre_N_k):length(post_N_k)-1) length(post_N_k)];
    end
    %elseif length(comps) == 2 %merge
    %    pre_posterior.gamma(:,1) = pre_gamma(:,pre_idx==comps(1));
    %    pre_posterior.gamma(:,2) = pre_gamma(:,pre_idx==comps(2));
    %    if comps(1) > comps(2)
    %        merged_comp = comps(1) - 1;
    %    else
    %        merged_comp = comps(1);
    %    end
    %    post_posterior.gamma(:,1) = post_gamma(:,post_idx==merged_comp);
    %end
    
    
    
    pre_S_sk = eval_S_sk(subdata,pre_posterior,prior,options,pre_comps);
    post_S_sk = eval_S_sk(subdata,post_posterior,prior,options,post_comps);
    
    lse_pre_S_sk = helper_logsumexp(pre_S_sk,max(pre_S_sk,[],2));
    lse_post_S_sk = helper_logsumexp(post_S_sk,max(post_S_sk,[],2));
    
    num_singlet = size(subdata.singlets,2);
    
    delta_fe = -sum(post_component_KL) + sum(pre_component_KL) ...
               -sum(post_stick_KL) + sum(pre_stick_KL) ...
               +mag_factor*sum(lse_post_S_sk(1:num_singlet)) ...
               -mag_factor*sum(lse_pre_S_sk(1:num_singlet));
           
    if ~isempty(subdata.Nc)
        delta_fe = delta_fe ...
                 + mag_factor*subdata.Nc*lse_post_S_sk(num_singlet+1:end) ...
                 - mag_factor*subdata.Nc*lse_pre_S_sk(num_singlet+1:end);
    end
    
end                              