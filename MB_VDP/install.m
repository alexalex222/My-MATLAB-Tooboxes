%% Installation script (just compiles the MEX functions)
mex -O helper_B_diag.c
mex -O helper_B_full.c
mex -O helper_ELPX.c
mex -O helper_ELPX_clump.c
mex -O helper_ELPX_full.c
mex -O helper_column_sub.c
mex -O helper_logsumexp.c
