NUM_CASES = 139;

h = waitbar(0,'Please wait...');
stats_per_case = zeros(NUM_CASES, 2);
for case_number=1:NUM_CASES
    case_signal = all_data{case_number, 1};
    all_parameter_hat = parameters_per_case{case_number};
    num_voxels = size(case_signal, 1);
    noise_variances = zeros(num_voxels, 1);
    for j=1:num_voxels
        noise_variances(j) = evar(case_signal(j, :));
    end
    sorted_desc_nv = sort(noise_variances, 'descend');
    clip_five = round(size(sorted_desc_nv, 1)*0.05);
    [slope, intercept, MSE, R2, S] = logfit(1:clip_five, sorted_desc_nv(1:clip_five), 'loglog');
    stats_per_case(case_number, 1) = R2;
    clip_one = round(size(sorted_desc_nv, 1)*0.01);
    [slope, intercept, MSE, R2, S] = logfit(1:clip_one, sorted_desc_nv(1:clip_one), 'loglog');
    stats_per_case(case_number, 2) = R2;
    waitbar(case_number / NUM_CASES);
end
close(h);