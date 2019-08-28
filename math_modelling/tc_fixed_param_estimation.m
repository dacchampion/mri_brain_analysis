% Set options for the non­linear fitting algorithm:
h=optimset('MaxFunEvals', 20000,'Algorithm','levenberg-marquardt','TolX',1e-6, 'Display', 'off');

return_data = LLS_PV(all_data);
num_cases = 1;

S0_per_case = zeros(num_cases, 1);
m1_per_case = zeros(num_cases, 1);
m2_per_case = zeros(num_cases, 1);
m3_per_case = zeros(num_cases, 1);
RESNORMS = zeros(num_cases, 1);
for case_number=1:num_cases
    t_start = tic;
    S0_vec = return_data{case_number, 1}';
    
    TE_values = all_data{case_number, 2};
    
    case_voxels = all_data{case_number, 1};
    startx_mat = [sqrt(S0_vec), ones(size(S0_vec, 1), 1)*0.13, ones(size(S0_vec, 1), 1)*0.42, ones(size(S0_vec, 1), 1)*0.35];
    num_vox = size(case_voxels, 1);
    all_parameter_hat = zeros(num_vox, 4);
    average_RESNORM = 0;
    for vox_ind=1:num_vox
        true_S = case_voxels(vox_ind, :);
        startx = startx_mat(vox_ind, :);

        [parameter_hat, RESNORM] = fminunc('three_compartment', startx, h, true_S, TE_values);
        average_RESNORM = average_RESNORM + RESNORM;
        all_parameter_hat(vox_ind, 1) = parameter_hat(1)^2;
        all_parameter_hat(vox_ind, 2) = sin(parameter_hat(2))^2;
        all_parameter_hat(vox_ind, 3) = (1-all_parameter_hat(vox_ind, 2))*sin(parameter_hat(3))^2;
        all_parameter_hat(vox_ind, 4) = 1-all_parameter_hat(vox_ind, 2)-all_parameter_hat(vox_ind, 3);
    end
    RESNORMS(case_number) = average_RESNORM/num_vox;
    S0_per_case(case_number) = mean(all_parameter_hat(:, 1));
    m1_per_case(case_number) = mean(all_parameter_hat(:, 2));
    m2_per_case(case_number) = mean(all_parameter_hat(:, 3));
    m3_per_case(case_number) = mean(all_parameter_hat(:, 4));
    t_elapsed = toc(t_start);
    fprintf(strcat('Time elapsed: ', num2str(t_elapsed), ', for case', num2str(case_number),'\n'));
end