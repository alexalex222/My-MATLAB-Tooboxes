function [free_energy, posterior, q_z] = iterate_posterior(data,posterior,prior,options,iter,do_sort,suppress_output)
    
    if nargin < 6;
        do_sort = true;
        iter = inf;
        suppress_output = false;
    end
    
    if ~suppress_output
        disp('Updating until convergence...');
    end
    
    free_energy = -inf;
    i = 0;
    start_sort = 0;
    
    while true
        i = i+1;
        [new_free_energy, S_sk] = eval_free_energy(data, posterior, prior, options);
        if ~suppress_output
            disp_status(new_free_energy, posterior);
        end
        % if the free energy isn't improved
        if (isinf(iter) && ~free_energy_improved(free_energy, new_free_energy, options) ) ...
           || (~isinf(iter) && i >=iter)
            free_energy = new_free_energy;
            if ~start_sort && do_sort %begin sorting process
                start_sort = 1;
            else %if post-sorting updating has completed
                break;
            end
        end

        free_energy = new_free_energy;
        q_z = eval_q_z(data, posterior, prior, options, S_sk);

        % check if the infinite part is small enough
        if sum(q_z.singlets(:,end))+sum(q_z.clumps((data.Nc>0),end)) > 1.0e-20
            q_z.singlets(:,end+1) = 0;  % if not make a new component
            q_z.clumps(:,end+1)   = 0;
        end
        if start_sort
            disp('Sorting...');
            N_k = sum(q_z.singlets,1) + data.Nc*q_z.clumps;
            [foo,I] = sort(N_k,2,'descend');
            q_z.singlets = q_z.singlets(:,I);
            q_z.clumps   = q_z.clumps(:,I);
        end
        
        % if last component is very small
        if (sum(q_z.singlets(:,end-1))+sum(q_z.clumps((data.Nc>0),end-1))) < 1.0e-10
            q_z.singlets(:,end-1) = []; %get rid of it
            q_z.clumps(:,end-1)   = [];
        end
        
        posterior = eval_posterior(data, q_z, prior, options);
    end
    if suppress_output
        disp_status(new_free_energy, posterior);
    end
    
end
%%
function disp_status(free_energy, posterior)

    N_k = posterior.true_N_k;
    disp(['F=' num2str(free_energy) '; K= ' num2str(length(N_k)) '; N_k=[' num2str(N_k, ' %0.5g ') '];'])

end