folder_name = '~/raw-data/cmbi_data/';
load(strcat(folder_name, 'info.mat'));

% Retrieve T2 file names
fileList = dir(strcat(folder_name, '*qt2.nii.gz'));
names_qt2_reg = {fileList.name};

% Retrieve mask file names
fileList = dir(strcat(folder_name, '*mask.nii.gz'));
names_masks = {fileList.name};

% Retrieve segmentation file names
fileList = dir(strcat(folder_name, '*seg.nii.gz'));
names_seg = {fileList.name};

% Retrieve parcellation file names
fileList = dir(strcat(folder_name, '*par1.nii.gz'));
names_par1 = {fileList.name};

% Retrieve simplified parcellation file names
fileList = dir(strcat(folder_name, '*par2.nii.gz'));
names_par2 = {fileList.name};

% load the echo times
TEs = [13, 16, 20, 25, 30, 40, 50, 85, 100, 150];

% store data about the images (with the TEs) in a cell array
% we have 5 files per image (ignoring the mask)
% we have numel(names_qt2_reg) = 6 cases
all_data = cell(numel(names_qt2_reg), 6);

for num = 1:numel(names_qt2_reg)
    t_start = tic;
    % collect file directories
    file_dir0 = strcat(folder_name, names_masks{num});
    file_dir2 = strcat(folder_name, names_qt2_reg{num});
    file_dir3 = strcat(folder_name, names_seg{num});
    file_dir4 = strcat(folder_name, names_par1{num});
    file_dir5 = strcat(folder_name, names_par2{num});

    
    % load mask
    try
        data = load_nii(file_dir0);
    catch
        data = load_untouch_nii(file_dir0);
    end
    mask = data.img;
    mask_vec = mask(:);
    
    % load the image
    try
        data = load_nii(file_dir2);
    catch
        data = load_untouch_nii(file_dir2);
    end
    
    qt2_reg_RAW = data.img; % size (96, 96, 55, num_TEs)
    num_TEs = size(qt2_reg_RAW, 4);
    qt2_reg_RAW_matrix = reshape(qt2_reg_RAW, 96*96*51, num_TEs);

    % load the segmentation
    try
        data = load_nii(file_dir3);
    catch
        data= load_untouch_nii(file_dir3);
    end
    seg_RAW = data.img; % size (96, 96, 55, 6)
    seg_RAW_matrix = reshape(seg_RAW, 96*96*51, 6);

    % load par1 (parcellation)
    try
        data = load_nii(file_dir4);
    catch
        data = load_untouch_nii(file_dir4);
    end
    par1_RAW= data.img; % size (96, 96, 55)
    par1_RAW_matrix = reshape(par1_RAW, 96*96*51, 1);

    % load par2 (simplified parcellation)
    try
        data = load_nii(file_dir5);
    catch
        data = load_untouch_nii(file_dir5);
    end
    par2_RAW= data.img; % size (96, 96, 55)
    par2_RAW_matrix = reshape(par2_RAW, 96*96*51, 1);
    
    % apply mask to data
    num_voxels_b4_mask = size(qt2_reg_RAW_matrix, 1);
    for itm=1:num_voxels_b4_mask
        if mask_vec(itm)~=0 && sum(qt2_reg_RAW_matrix(itm,:)<0)>0
            mask_vec(itm) = 0;
        end
    end    
    qt2_reg = zeros(nnz(mask_vec), num_TEs);
    seg = zeros(nnz(mask_vec), 6);
    par1 = zeros(nnz(mask_vec), 1);
    par2 = zeros(nnz(mask_vec), 1);
    
    
    row=0;
    for itm=1:num_voxels_b4_mask
        if mask_vec(itm)~=0
            row = row + 1;
            qt2_reg(row, :) = qt2_reg_RAW_matrix(itm, :);
            segmentation_signal = seg_RAW_matrix(itm, :);
            segmentation_signal(segmentation_signal<0)=0;
            segmentation_signal(segmentation_signal>1)=1;
            seg(row, :) = segmentation_signal;
            par1(row, :) = par1_RAW_matrix(itm, :);
            par2(row, :) = par2_RAW_matrix(itm, :);
        end
    end
    
    % note that the mask gave us has nnz(mask_vec) non-zero elements
    all_data{num, 1} = qt2_reg; % store the resulting image of size (nnz(mask_vec), num_TEs)
    all_data{num, 2} = seg;
    all_data{num, 3} = par1;
    all_data{num, 4} = par2;
    all_data{num, 5} = reshape(mask_vec, [size(mask, 1), size(mask, 2), size(mask, 3)]);
    all_data{num, 6} = [cohort(num, :), stats(num, 2:4)];
    t_elapsed = toc(t_start);
    fprintf(strcat('Time elapsed: ', num2str(t_elapsed), ', for case', num2str(num),'\n'));
end