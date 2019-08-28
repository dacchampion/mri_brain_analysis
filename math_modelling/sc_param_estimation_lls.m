NUM_CASES = 1;

lls_params = LLS_PV(all_data);
RESNORM = zeros(NUM_CASES, 1);
for case_number=1:NUM_CASES
    t_start = tic;
    A = all_data{case_number, 1};
    TEs = all_data{case_number, 2};
    S0 = lls_params{case_number, 1};
    T2 = lls_params{case_number, 2};    
    total_vox = size(A, 1);
    for vox_num=1:total_vox
        S = S0(vox_num)*exp(-TEs/T2(vox_num));
        RESNORM(case_number) = RESNORM(case_number) + sum((A(vox_num, :) - S).^2);
    end
    RESNORM(case_number) = RESNORM(case_number)/total_vox;
    t_elapsed = toc(t_start);
    fprintf(strcat('Time elapsed: ', num2str(t_elapsed), ', for case ', num2str(case_number),'\n'));
end