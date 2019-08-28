case_number = 116;
case_signal = all_data{case_number, 1};

%TEs = all_data{case_number, 2};
all_parameter_hat = parameters_per_case{case_number};

num_voxels = size(case_signal, 1);
noise_variances = zeros(num_voxels, 1);
for i=1:num_voxels
    noise_variances(i) = evar(case_signal(i, :));
end

[sorted_nv, unsorted_idxs] = sort(noise_variances);

mid_voxel = round(num_voxels/2);
high_var = sorted_nv(num_voxels);
mid_var = sorted_nv(mid_voxel);
low_var = sorted_nv(1);

high_noise_s = case_signal(unsorted_idxs(num_voxels), :);
mid_noise_s = case_signal(unsorted_idxs(mid_voxel), :);
low_noise_s = case_signal(unsorted_idxs(1), :);

hv_params = all_parameter_hat(unsorted_idxs(num_voxels), :);
mv_params = all_parameter_hat(unsorted_idxs(mid_voxel), :);
lv_params = all_parameter_hat(unsorted_idxs(1), :);

[hv_resnorm, estimated_hv_s] = three_compartement_untrans(hv_params, high_noise_s, TEs);
[mv_resnorm, estimated_mv_s] = three_compartement_untrans(mv_params, mid_noise_s, TEs);
[lv_resnorm, estimated_lv_s] = three_compartement_untrans(lv_params, low_noise_s, TEs);

noisy_signal = high_noise_s + sqrt(high_var).*randn(1, size(TEs, 2));
hv_snr = snr(estimated_hv_s, high_noise_s-noisy_signal);
noisy_signal = mid_noise_s + sqrt(mid_var).*randn(1, size(TEs, 2));
mv_snr = snr(estimated_mv_s, mid_noise_s-noisy_signal);
noisy_signal = low_noise_s + sqrt(low_var).*randn(1, size(TEs, 2));
lv_snr = snr(estimated_lv_s, low_noise_s-noisy_signal);

idx_for_text = size(TEs, 2);
subplot(1,3,1), plot(TEs, low_noise_s, '-o', TEs, estimated_lv_s, '-x'), title('Estimation with low-variance noise')
snr_text =strcat('SNR \rightarrow ', num2str(lv_snr));
sigma_text = strcat('Noise \sigma \rightarrow ', num2str(sqrt(low_var)));
text(TEs(1, 6), estimated_lv_s(1, idx_for_text-4), snr_text);
text(TEs(1, 5), estimated_lv_s(1, idx_for_text-5), sigma_text);
legend('Signal', 'Estimation')
xlabel('Echo Time');
ylabel('Intensity');
subplot(1,3,2), plot(TEs, mid_noise_s, '-o', TEs, estimated_mv_s, '-x'), title('Estimation with mid-variance noise')
snr_text =strcat('SNR \rightarrow ', num2str(mv_snr));
sigma_text = strcat('Noise \sigma \rightarrow ', num2str(sqrt(mid_var)));
text(TEs(1, 8), estimated_mv_s(1, idx_for_text-2), snr_text);
text(TEs(1, 7), estimated_mv_s(1, idx_for_text-3), sigma_text);
legend('Signal', 'Estimation')
xlabel('Echo Time');
ylabel('Intensity');
subplot(1,3,3), plot(TEs, high_noise_s, '-o', TEs, estimated_hv_s, '-x'), title('Estimation with high-variance noise')
snr_text = strcat('SNR \rightarrow ', num2str(hv_snr));
sigma_text = strcat('Noise \sigma \rightarrow ', num2str(sqrt(high_var)));
text(TEs(1, 9), estimated_hv_s(1, idx_for_text-1), snr_text);
text(TEs(1, 6), estimated_hv_s(1, idx_for_text-4), sigma_text);
legend('Signal', 'Estimation')
xlabel('Echo Time');
ylabel('Intensity');