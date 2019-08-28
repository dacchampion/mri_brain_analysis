function [return_data] = LLS_PV_v2(all_data, TEs_raw)
% Linear least squares function
% all data is a cell array


num_cases = size(all_data, 1);
return_data = cell(num_cases, 2);

for img=1:num_cases
    y_mat = cell2mat(all_data(img, 1));
    y_mat(y_mat<0) = 1; % All negative values to 1 TEMP!!!!!!!!!!!!!!
    y_mat = log(y_mat);
    
    % Get A by repeating TEs in a column vector
    num_TEs = size(TEs_raw, 1);

    A = [ones(num_TEs, 1) TEs_raw];
    b_part = pinv(A' * A) * (A');
    b = b_part * (y_mat'); % size = (2, num_voxels)
    
    S0 = exp(b(1, :));
    T2 = -1./(b(2, :));
    
    return_data{img, 1} = S0;
    return_data{img, 2} = T2;
end