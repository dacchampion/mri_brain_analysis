function [return_data] = reconstruct_image_sequence(image_obs_per_TE, image_mask)
% When reconstructing images, we set all the values not included in the
% mask to zero, and only place voxels in the positions for which the mask
% was 1
% INPUTS::
% image_obs_per_TE = a cell array size (num_cases, 1), where each element has size (num_masked_voxels, num_TEs)
% Note: num_masked_voxels & num_TEs are different for each element of the cell array
% image_mask = the mask for each case

% useful variable
num_cases = size(image_mask, 1);

% store data to be returned in this cell array
return_data = cell(num_cases, 1);


for img=1:num_cases
    % corresponds to the to the voxel from the masked input image
    next_vox = 1;

    % size (num_masked_voxels, num_TE)
    input_image = image_obs_per_TE{img};
    num_TE = size(input_image, 2); % size of 4th dimension

    % size (96*96*55, 1)
    input_mask = image_mask{img};
    mask1 = input_mask(:);
    input_mask_vec = repmat(mask1, num_TE, 1); % a repeated mask for every TE

    % total voxel in reconstructed image
    total_vox = size(input_mask_vec, 1);
    reconstructed = zeros(total_vox, 1);

    for ii=1:total_vox

        if input_mask_vec(ii) == 1
            reconstructed(ii) = input_image(next_vox);
            next_vox = next_vox +1;
        end
    end
    return_data{img} = reshape(reconstructed, [size(input_mask, 1), size(input_mask, 2), size(input_mask, 3), num_TE]);
    
end