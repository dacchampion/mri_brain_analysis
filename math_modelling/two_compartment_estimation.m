%% multi compartment with constraints

first_estimation = LLS_PV(all_data);

% Set options for the non­linear fitting algorithm:
h=optimset('MaxFunEvals', 20000,'Algorithm','levenberg-marquardt','TolX',1e-6, 'Display', 'off');

num_cases = 1;

% startx_mat = [mean(sqrt(S0_vec)) * ones(size(S0_vec)), sqrt(T2s), m1];

S0_per_case = zeros(num_cases, 1);
T2_1_per_case = zeros(num_cases, 1);
T2_2_per_case = zeros(num_cases, 1);
m1_per_case = zeros(num_cases, 1);
RESNORMS = zeros(num_cases, 1);
for case_number=1:num_cases
    S0_vec = first_estimation{case_number, 1}';    
    TE_values = all_data{case_number, 2};
    true_S_mat = all_data{case_number, 1};
    num_vox = size(true_S_mat, 1);
    m1 = 0.5*ones(num_vox, 1);
    T2s = [20*ones(num_vox, 1), 200*ones(num_vox, 1)];
    startx_mat = [sqrt(S0_vec), sqrt(T2s), m1];
    all_parameter_hat = zeros(num_vox, 3);
    average_RESNORM = 0;
    t_start = tic;
    for vox_ind=1:num_vox
        true_S = true_S_mat(vox_ind, :);
        startx = startx_mat(vox_ind, :);

        [parameter_hat, RESNORM] = fminunc('two_compartment',startx, h, true_S, TE_values);
        average_RESNORM = average_RESNORM + RESNORM;
        all_parameter_hat(vox_ind, 1:3) = parameter_hat(1:3).^2;
        all_parameter_hat(vox_ind, 4) = (sin(parameter_hat(4)))^2;
    end
    RESNORMS(case_number) = average_RESNORM/num_vox;
    S0_per_case(case_number) = mean(all_parameter_hat(:, 1));
    T2_1_per_case(case_number) = mean(all_parameter_hat(:, 2));
    T2_2_per_case(case_number) = mean(all_parameter_hat(:, 3));
    m1_per_case(case_number) = mean(all_parameter_hat(:, 4));    
    t_elapsed = toc(t_start);
    fprintf(strcat('Time elapsed: ', num2str(t_elapsed), ', for case ', num2str(case_number),' with RESNORM=', num2str(RESNORMS(case_number)), '\n'));
end


% Initialise: m1=m2=0.5; T2 = [80, 200]; S0 from non-linear est section
% e = 902.6600, average_RESNORM = 5.0559e+04 i.e. slower but has lower RESNORM

% Using m1=m2=0.5; T2 = [80, 200]; S0 set to mean of non-linear est section
% e = X, average_RESNORM = Y i.e. Z


% Using m1=GM+WM segments; T2 = [80, 200]; S0 set to mean of non-linear est section
% e = X, average_RESNORM = Y i.e. Z

% Using m1=CSF segments; T2 = [80, 200]; S0 set to mean of non-linear est section
% e = X, average_RESNORM = Y i.e. Z