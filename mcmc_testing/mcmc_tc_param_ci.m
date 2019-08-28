MCMC_SAMPLES = 100000;
BURNIN = 50000;

EP_CASES = [20, 137, 70, 54, 81, 97, 84, 10, 124, 107, 111, 21, 100, 73, 96, 76, 17, 79, 11, 89];
FP_CASES = [14, 15, 16, 18, 32, 36, 37, 42, 47, 52, 56, 58, 69, 86, 92, 94, 98, 99, 102, 105];

mean_S0 = 0;
mean_m1 = 0;
mean_m2 = 0;
num_cases4_analysis = size(FP_CASES, 2);
ep_data = cell(num_cases4_analysis, 1);
for idx=1:num_cases4_analysis
    case_number = FP_CASES(idx);
    case_parameter = parameters_per_case{case_number};
    mean_S0 = mean_S0 + mean(case_parameter(:, 1));
    mean_m1 = mean_m1 + mean(case_parameter(:, 2));
    mean_m2 = mean_m2 + mean(case_parameter(:, 3));
    
    case_voxels = all_data{case_number, 1};
    segments = all_data{case_number, 2};
    wm_segment = segments(:, 4);
    
    wm_voxels = zeros(nnz(wm_segment), size(TEs, 2));
    ii = 1;
    for j=1:size(case_voxels, 1)
        if wm_segment(j)~=0
            wm_voxels(ii, :) = case_voxels(j, :);
            ii = ii +1;
        end
    end
    ep_data{idx} = wm_voxels;
end

mean_S0 = mean_S0/num_cases4_analysis;
mean_m1 = mean_m1/num_cases4_analysis;
mean_m2 = mean_m2/num_cases4_analysis;

A_hat = [];
for k=1:num_cases4_analysis
    A_hat = vertcat(A_hat, ep_data{k});
end
data{1} = TEs;
data{2} = A_hat;

likelihood = @logLikelihood;
model = @three_compartment_model;

%prior = {'S0', 'gaussian', mean_S0, 10, ''; ...
%         'm1', 'gaussian', mean_m1, 0.05, ''; ...
%         'm2', 'gaussian', mean_m2, 0.1, ''};

prior = {'S0', 'gaussian', mean_S0, mean_S0*0.2, ''; ...
         'm1', 'uniform', mean_m1, mean_m1*0.1, 'cyclic'; ...
         'm2', 'uniform', mean_m2, mean_m2*0.1, 'cyclic'};    
extraparams = {'m3', 1-mean_m1-mean_m2};

[x_t, logP] = mcmc_sampler(data, ...
                           likelihood, ... 
                           model, ...
                           prior, ...
                           extraparams, ...
                           'recalcprop', 5000, ...
                           'Nmcmc', MCMC_SAMPLES, ...
                           'Nburnin', BURNIN);