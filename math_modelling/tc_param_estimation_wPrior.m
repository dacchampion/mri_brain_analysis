% Set options for the non­linear fitting algorithm:
h=optimset('MaxFunEvals', 20000,'Algorithm','levenberg-marquardt','TolX',1e-6, 'Display', 'off');

first_estimation = LLS_PV(all_data);

%first_estimation = LLS_PV_v2(all_data, TEs');
num_cases = 6;

S0_per_case = zeros(num_cases, 1);
m1_per_case = zeros(num_cases, 1);
m2_per_case = zeros(num_cases, 1);
parameters_per_case = cell(num_cases, 1);
RESNORMS = zeros(num_cases, 1);
for case_number=6:num_cases
    t_start = tic;
    S0_vec = first_estimation{case_number, 1}';
    
    case_voxels = all_data{case_number, 1};
    TE_values = all_data{case_number, 2};
    %TE_values = TEs;
    segments = all_data{case_number, 3};
    %segments = all_data{case_number, 2};

    num_vox = size(case_voxels, 1);
    all_parameter_hat = zeros(num_vox, 3);
    average_RESNORM = 0;
    for vox_ind=1:num_vox
        true_S = case_voxels(vox_ind, :);
        
        s0 = sqrt(S0_vec(vox_ind));
        
        gm = segments(vox_ind, 3);
        wm = segments(vox_ind, 4);
        bg = segments(vox_ind, 5);
        bs = segments(vox_ind, 6);
        
        mwf = wm*0.25 +gm*0.03 + bg*0.01 + bs*0.01;
        tf = wm*0.55 +gm*0.3 + bg*0.1 + bs*0.1;
        csf = 1-mwf-tf;
        startx = [s0, mwf, tf, csf];
        try
            [parameter_hat, RESNORM] = fminunc('three_compartment', startx, h, true_S, TE_values);
        catch ME1
            mwf = rand;
            tf = sin(mwf)^2;
            csf = 1-mwf-tf;
            startx = [s0, mwf, tf, csf];
            try
                [parameter_hat, RESNORM] = fminunc('three_compartment', startx, h, true_S, TE_values);
            catch ME2
                  fprintf('In SECOND catch, case number = %d', case_number);
                  rethrow(ME1)
            end
        end
            
        average_RESNORM = average_RESNORM + RESNORM;
        all_parameter_hat(vox_ind, 1) = parameter_hat(1)^2;
        all_parameter_hat(vox_ind, 2) = sin(parameter_hat(2))^2;
        all_parameter_hat(vox_ind, 3) = (1-all_parameter_hat(vox_ind, 2))*sin(parameter_hat(3))^2;
    end
    RESNORMS(case_number) = average_RESNORM/num_vox;
    S0_per_case(case_number) = mean(all_parameter_hat(:, 1));
    m1_per_case(case_number) = mean(all_parameter_hat(:, 2));
    m2_per_case(case_number) = mean(all_parameter_hat(:, 3));
    parameters_per_case{case_number}= all_parameter_hat;
    t_elapsed = toc(t_start);
    fprintf(strcat('Time elapsed: ', num2str(t_elapsed), ', for case ', num2str(case_number),' with RESNORM=', num2str(RESNORMS(case_number)), '\n'));
    clear('case_voxels', 'average_RESNORM', 'bg', 'bs', 'csf', 'gm', 'mwf', 'num_cases', 'num_vox', 'parameter_hat', 'RESNORM', 's0', 'S0_vec', 'segments', 'startx', 't_elapsed', 't_start', 'TE_values', 'tf', 'true_S', 'vox_ind', 'wm');
    save('checkpoint.mat', 'all_parameter_hat', 'case_number', 'S0_per_case', 'm1_per_case', 'm2_per_case', 'RESNORMS', 'parameters_per_case');
end