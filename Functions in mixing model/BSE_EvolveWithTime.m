function[N_parent_BSE,N_daughter_BSE,N_stable_BSE,R_BSE] = ...
    BSE_EvolveWithTime(t,nt,N_parent_BSE_t0,N_daughter_BSE_t0,...
    N_stable_BSE_t0,R_BSE_t0,lambda_parent)

N_parent_BSE = nan(size(t)); N_parent_BSE(1) = N_parent_BSE_t0; 
N_daughter_BSE = nan(size(t)); N_daughter_BSE(1) = N_daughter_BSE_t0;
N_stable_BSE = nan(size(t)); N_stable_BSE(1) = N_stable_BSE_t0;
R_BSE = nan(size(t)); R_BSE(1) = R_BSE_t0;
for i = 2:nt
    N_parent_BSE(i) = N_parent_BSE(1) * exp(-lambda_parent*t(i));
    N_daughter_BSE(i) = N_daughter_BSE(1) + N_parent_BSE(1) * (1 - exp(-lambda_parent*t(i)));
    N_stable_BSE(i) = N_stable_BSE(1);
    R_BSE(i) = N_daughter_BSE(i) / N_stable_BSE(i);
end