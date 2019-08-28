EP_CASES = [20, 137, 70, 54, 81, 97, 84, 10, 124, 107, 111, 21, 100, 73, 96, 76, 17, 79, 11, 89];
FT_CASES = [14, 15, 16, 18, 32, 36, 37, 42, 47, 52, 56, 58, 69, 86, 92, 94, 98, 99, 102, 105];

NUM_PARCELS = 12;
num_cases = size(EP_CASES, 2);
ft_mwf_per_parcel = cell(NUM_PARCELS, 1);
ep_mwf_per_parcel = cell(NUM_PARCELS, 1);
ft_tf_per_parcel = cell(NUM_PARCELS, 1);
ep_tf_per_parcel = cell(NUM_PARCELS, 1);
for i=1:num_cases
    case_number = FT_CASES(i);
    case_parameters = parameters_per_case{case_number};
    mwf = case_parameters(:, 2);
    case_parcellation1 = all_data{case_number, 3};
    num_voxels = size(mwf, 1);
    for j=1:NUM_PARCELS
        labeled_mwf = mwf(case_parcellation1==j);
        parcel_vector = ft_mwf_per_parcel{j};
        parcel_vector = vertcat(parcel_vector, labeled_mwf);
        ft_mwf_per_parcel{j} = parcel_vector;
    end
    
    case_number = EP_CASES(i);
    case_parameters = parameters_per_case{case_number};
    mwf = case_parameters(:, 2);
    case_parcellation1 = all_data{case_number, 3};
    num_voxels = size(mwf, 1);    
    for j=1:NUM_PARCELS
        labeled_mwf = mwf(case_parcellation1==j);
        parcel_vector = ep_mwf_per_parcel{j};
        parcel_vector = vertcat(parcel_vector, labeled_mwf);
        ep_mwf_per_parcel{j} = parcel_vector;
    end
    
    case_number = FT_CASES(i);
    case_parameters = parameters_per_case{case_number};
    tf = case_parameters(:, 3);
    case_parcellation1 = all_data{case_number, 3};
    num_voxels = size(tf, 1);    
    for j=1:NUM_PARCELS
        labeled_tf = tf(case_parcellation1==j);
        parcel_vector = ft_tf_per_parcel{j};
        parcel_vector = vertcat(parcel_vector, labeled_tf);
        ft_tf_per_parcel{j} = parcel_vector;
    end
    
    case_number = EP_CASES(i);
    case_parameters = parameters_per_case{case_number};
    tf = case_parameters(:, 3);
    case_parcellation1 = all_data{case_number, 3};
    num_voxels = size(tf, 1);
    for j=1:NUM_PARCELS
        labeled_tf = tf(case_parcellation1==j);
        parcel_vector = ep_tf_per_parcel{j};
        parcel_vector = vertcat(parcel_vector, labeled_tf);
        ep_tf_per_parcel{j} = parcel_vector;
    end    
end

% Post processing each parcel to contain the same number of observations
% for i=1:NUM_PARCELS
%     ep_parcel = ep_mwf_per_parcel{i};
%     ft_parcel = ft_mwf_per_parcel{i};
%     clipping_idx = min(size(ep_parcel, 1), size(ft_parcel, 1));
%     ep_mwf_per_parcel{i} = ep_parcel(1:clipping_idx, 1);
%     ft_mwf_per_parcel{i} = ft_parcel(1:clipping_idx, 1);
% end