function q_z = eval_q_z(data, posterior, prior, options, S_sk)
    num_singlet = size(data.singlets,2);
    
    if nargin < 5
        S_sk = eval_S_sk(data, posterior, prior, options);
    end
    
    lse_S_sk = helper_logsumexp(S_sk,max(S_sk,[],2));
    helper_column_sub(S_sk,lse_S_sk);
    
    q_z.singlets = exp(S_sk(1:num_singlet,:));
    q_z.clumps   = exp(S_sk(num_singlet+1:end,:));

end