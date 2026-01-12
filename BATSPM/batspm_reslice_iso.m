function batspm_reslice_iso(nii_path, target_res, force_z_thick)
% BATSPM_RESLICE_ISO (AUTO-REPAIR + SAFETY CLAMP)
% 1. Checks for "Swapped Header" (Pencil Effect) and fixes it.
% 2. Finds the Brain (Anchor).
% 3. Clamps huge ghosts (>50mm).
% 4. Reslices to isotropic.

    if nargin < 2, target_res = 0.4; end 

    % 1. Load Data
    V = spm_vol(nii_path);
    Y = double(spm_read_vols(V)); 
    
    fprintf('Processing %s...\n', spm_file(nii_path, 'filename'));
    
    % --- STEP 1.5: HEADER AUTO-REPAIR (The "Pencil" Fix) ---
    % We check if the smallest dimension (slices) has a tiny voxel size.
    % If so, we swap it with the largest voxel size.
    
    dims = V.dim;
    % Calculate current voxel sizes from the matrix columns
    vox_sizes = sqrt(sum(V.mat(1:3,1:3).^2));
    
    [min_dim_val, min_dim_idx] = min(dims);       % e.g., 18 (Index 1)
    [max_vox_val, max_vox_idx] = max(vox_sizes);  % e.g., 0.8 (Index 3)
    
    % If the fewest slices (18) don't have the thickest size (0.8)...
    if min_dim_idx ~= max_vox_idx
        fprintf('   -> DETECTED HEADER MISMATCH (Pencil Effect).\n');
        fprintf('      Swapping voxel sizes: Axis %d (was %.2f) <-> Axis %d (was %.2f)\n', ...
            min_dim_idx, vox_sizes(min_dim_idx), max_vox_idx, vox_sizes(max_vox_idx));
        
        % REPAIR STRATEGY:
        % 1. Extract rotation vectors (Direction)
        R = V.mat(1:3, 1:3);
        norm_R = R ./ vox_sizes; % Unit vectors
        
        % 2. Swap the scaling factors
        new_vox_sizes = vox_sizes;
        new_vox_sizes(min_dim_idx) = vox_sizes(max_vox_idx); % Assign 0.8 to Axis 1
        new_vox_sizes(max_vox_idx) = vox_sizes(min_dim_idx); % Assign 0.16 to Axis 3
        
        % 3. Rebuild the matrix
        V.mat(1:3, 1:3) = norm_R .* new_vox_sizes;
        
        fprintf('      New Physical FOV: %.1f x %.1f x %.1f mm\n', ...
            dims .* new_vox_sizes);
    end

    % --- OVERRIDE (Manual force if needed) ---
    if nargin > 2 && ~isempty(force_z_thick)
        V.mat(3,3) = force_z_thick * sign(V.mat(3,3)); 
    end

    % 2. FIND THE BRAIN (Anchor Method)
    Y_smooth = smooth3(Y, 'box', 3);
    [max_val, max_idx] = max(Y_smooth(:));
    mask = Y_smooth > (max_val * 0.15);
    
    CC = bwconncomp(mask);
    found_brain = false;
    for i = 1:CC.NumObjects
        if ismember(max_idx, CC.PixelIdxList{i})
            clean_mask = false(size(mask));
            clean_mask(CC.PixelIdxList{i}) = true;
            found_brain = true;
            break;
        end
    end
    if ~found_brain, clean_mask = mask; end
    
    % 3. CALCULATE BOUNDING BOX
    [ind_x, ind_y, ind_z] = ind2sub(size(clean_mask), find(clean_mask));
    pad = 5;
    min_ijk = [min(ind_x)-pad, min(ind_y)-pad, min(ind_z)-pad];
    max_ijk = [max(ind_x)+pad, max(ind_y)+pad, max(ind_z)+pad];
    
    % --- SAFETY CLAMP (50mm Limit) ---
    % Recalculate voxel size based on the FIXED header
    current_vox = sqrt(sum(V.mat(1:3,1:3).^2));
    box_size_vox = max_ijk - min_ijk;
    box_size_mm = box_size_vox .* current_vox;
    
    max_allowed_mm = 50; 
    
    for d = 1:3
        if box_size_mm(d) > max_allowed_mm
            % Clamp logic
            center_idx = round((min_ijk(d) + max_ijk(d)) / 2);
            half_width = round((max_allowed_mm / current_vox(d)) / 2);
            min_ijk(d) = center_idx - half_width;
            max_ijk(d) = center_idx + half_width;
        end
    end
    
    min_ijk = max(min_ijk, 1);
    max_ijk = min(max_ijk, [V.dim(1) V.dim(2) V.dim(3)]);
    
    % 4. GENERATE GRID & INTERPOLATE
    % (Standard Inverse Mapping)
    [crop_x, crop_y, crop_z] = ndgrid(min_ijk(1):max_ijk(1), ...
                                      min_ijk(2):max_ijk(2), ...
                                      min_ijk(3):max_ijk(3));
                                  
    corners_x = V.mat(1,1)*crop_x + V.mat(1,2)*crop_y + V.mat(1,3)*crop_z + V.mat(1,4);
    corners_y = V.mat(2,1)*crop_x + V.mat(2,2)*crop_y + V.mat(2,3)*crop_z + V.mat(2,4);
    corners_z = V.mat(3,1)*crop_x + V.mat(3,2)*crop_y + V.mat(3,3)*crop_z + V.mat(3,4);
    
    min_box = [min(corners_x(:)), min(corners_y(:)), min(corners_z(:))];
    max_box = [max(corners_x(:)), max(corners_y(:)), max(corners_z(:))];
    
    new_x = min_box(1) : target_res : max_box(1);
    new_y = min_box(2) : target_res : max_box(2);
    new_z = min_box(3) : target_res : max_box(3);
    
    [xq_mm, yq_mm, zq_mm] = ndgrid(new_x, new_y, new_z);
    
    F = griddedInterpolant(Y, 'linear', 'none');
    M_inv = inv(V.mat);
    
    src_x = M_inv(1,1)*xq_mm + M_inv(1,2)*yq_mm + M_inv(1,3)*zq_mm + M_inv(1,4);
    src_y = M_inv(2,1)*xq_mm + M_inv(2,2)*yq_mm + M_inv(2,3)*zq_mm + M_inv(2,4);
    src_z = M_inv(3,1)*xq_mm + M_inv(3,2)*yq_mm + M_inv(3,3)*zq_mm + M_inv(3,4);
    
    Y_new = F(src_x, src_y, src_z);
    Y_new(isnan(Y_new)) = 0;
    
    % 5. SAVE & CENTER
    V_new = V;
    V_new.dim = size(Y_new);
    V_new.mat = diag([target_res, target_res, target_res, 1]);
    
    total_mass = sum(Y_new(:));
    [gx, gy, gz] = ndgrid(1:V_new.dim(1), 1:V_new.dim(2), 1:V_new.dim(3));
    com_x = sum(Y_new(:) .* gx(:)) / total_mass;
    com_y = sum(Y_new(:) .* gy(:)) / total_mass;
    com_z = sum(Y_new(:) .* gz(:)) / total_mass;
    
    V_new.mat(1:3, 4) = -1 * (V_new.mat(1:3,1:3) * [com_x; com_y; com_z]);

    V_new.fname = nii_path; 
    spm_write_vol(V_new, Y_new);
    
    fprintf('   -> Success! Resliced to %.1f mm isotropic. (Dims: %d %d %d)\n', ...
        target_res, V_new.dim);
end