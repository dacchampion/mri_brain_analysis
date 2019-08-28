function [AICc, BIC] = model_criteria(N_per_voxel, tot_RESNORM, voxels, K)
    % K: number of all data points. K = all voxels*number of TEs
    % N: number of estimated parameters for one case
    N = (N_per_voxel+1)*voxels;
    AIC = 2*N + K*log(tot_RESNORM/K);
    AICc = AIC + (2*N*(N+1))/(K-N-1);
    BIC = N*log(K) + K*log(tot_RESNORM/K);
end