function image = reconstruct_image(input_image, image_mask)
    input_mask_vec = reshape(image_mask, [], 1);
    total_vox = size(input_mask_vec, 1);    
    reconstructed = zeros(total_vox, 1);
    next_vox = 1;
    for ii=1:total_vox
        if input_mask_vec(ii) == 1
            reconstructed(ii) = input_image(next_vox);
            next_vox = next_vox +1;
        end
    end
    image = reshape(reconstructed, [size(image_mask, 1), size(image_mask, 2), size(image_mask, 3)]);
end