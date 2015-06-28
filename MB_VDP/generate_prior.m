%Borrowed from Kenichi Kurihara's NIP 2006 implementation.
function hp_prior = generate_prior(data,alpha)
    hp_prior.xi0 = 0.01;
    eta_p = 1;
    
    if nargin > 1
        hp_prior.alpha = alpha;
    else
        hp_prior.alpha = 1;
    end
    
    [D,N] = size(data);
    
    covariance = cov(data');
    
    hp_prior.m0 = mean(data, 2);
    if D > 16
%        [dummy, max_eig] = power_method(covariance);
        max_eig = max(eig(covariance));
    else
        max_eig = max(eig(covariance));
    end
    hp_prior.eta0 = eta_p * D;
    hp_prior.B0 = hp_prior.eta0 * max_eig * eye(D) * hp_prior.xi0;
end