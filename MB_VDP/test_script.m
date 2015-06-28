%% 
load nist_features.mat
%% Full covariance
% The last argument is optional, it causes the program to output the
% current state to file after each learning round.  The string specifies
% the prefix of the files.
results = MB_VDP(nist_features,3000,3000,generate_prior(nist_features),'MB_3000_3000');
%% Diagonal covariance
prior = generate_prior(nist_features);
prior.B0 = diag(prior.B0);
results = MB_VDP(nist_features,3000,3000,prior,'MB_diag_3000_3000');