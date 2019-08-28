% Set options for the non­linear fitting algorithm:
h=optimset('MaxFunEvals', 20000,'Algorithm','levenberg-marquardt','TolX',1e-6, 'Display', 'off');

return_data = LLS_PV(all_data);
num_cases = 6;

S0_per_case = zeros(num_cases, 1);
m1_per_case = zeros(num_cases, 1);
RESNORMS = zeros(num_cases, 1);
for case_number=1:num_cases
    t_start = tic;
    S0_vec = return_data{case_number, 1}';
    T2_vec = return_data{case_number, 2}';
    
    TE_values = all_data{case_number, 2};
    
    case_voxels = all_data{case_number, 1};
    startx_mat = [sqrt(S0_vec), ones(size(S0_vec, 1), 1)*0.5];
    num_vox = size(case_voxels, 1);
    all_parameter_hat = zeros(num_vox, 4);
    average_RESNORM = 0;
    for vox_ind=1:num_vox
        true_S = case_voxels(vox_ind, :);
        startx = startx_mat(vox_ind, :);

        try
            [parameter_hat, RESNORM] = fminunc('two_constr_compartment', startx, h, true_S, TE_values);
        catch ME1
            S0 = startx(1)*rand;
            startx = [S0, 0.5];
            try
                [parameter_hat, RESNORM] = fminunc('two_constr_compartment', startx, h, true_S, TE_values);
            catch ME2
                fprintf('In SECOND catch, case number = %d', case_number);
                rethrow(ME1)
            end
        end            
        average_RESNORM = average_RESNORM + RESNORM;
        all_parameter_hat(vox_ind, 1) = parameter_hat(1)^2;
        all_parameter_hat(vox_ind, 2) = sin(parameter_hat(2))^2;
    end
    RESNORMS(case_number) = average_RESNORM/num_vox;
    S0_per_case(case_number) = mean(all_parameter_hat(:, 1));
    m1_per_case(case_number) = mean(all_parameter_hat(:, 4));
    t_elapsed = toc(t_start);
    fprintf(strcat('Time elapsed: ', num2str(t_elapsed), ', for case', num2str(case_number),'\n'));
end