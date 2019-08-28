function [return_data, RESNORMs] = nnLLS_PV_v2(all_data)
% Linear least squares function
% all_data is a cell array

num_cases = size(all_data, 1);
return_data = cell(num_cases, 2);
RESNORMs=zeros(num_cases,1);

for img=1:num_cases
    t_start = tic;
    RESNORM = 0;

    y_mat = all_data{img, 1}; % the all_data first column is the signal
    logY = log(y_mat);
    total_vox = size(logY, 1);
    
    % Get A by repeating TEs in a column vector
    TEs_raw = all_data{img, 2}';
    num_TEs = size(TEs_raw, 1);

    A = [ones(num_TEs, 1) TEs_raw];
    S0 = zeros(total_vox, 1);
    T2 = zeros(total_vox, 1);

    for vox_num=1:total_vox
        optnew = optimset('TolX',1e-4);
        b = lsqnonneg(A, logY(vox_num,:)', optnew);

        S0(vox_num) = exp(b(1));
        T2(vox_num) = -1/(b(2));

        % Compute the sum of square differences
        S = S0(vox_num)*exp(-TEs_raw/T2(vox_num));
        RESNORM = RESNORM + sum((y_mat(vox_num, :)' - S).^2);
    end
    return_data{img, 1} = S0;
    return_data{img, 2} = T2;
    RESNORMs(img) = RESNORM/total_vox;
    t_elapsed = toc(t_start);
    fprintf(strcat('Time elapsed: ', num2str(t_elapsed), ', for case ', num2str(img),'\n'));
end