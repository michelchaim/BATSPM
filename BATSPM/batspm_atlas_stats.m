function batspm_atlas_stats(jac_path, atlas_path, out_dir, bat_name)
% BATSPM_ATLAS_STATS (INTEGER-SAFE VERSION)
% 1. Reslices Atlas to Jacobian space using Nearest Neighbor (No decimals!).
% 2. Extracts stats for strictly integer regions (1, 2, 3...).

    fprintf('---------------------------------------------------\n');
    fprintf('Running Atlas Statistics for: %s\n', bat_name);

    if ~exist(jac_path, 'file'), error('Jacobian missing: %s', jac_path); end
    if ~exist(atlas_path, 'file'), error('Atlas missing: %s', atlas_path); end

    % --- STEP 1: RESLICE ATLAS (STRICT NEAREST NEIGHBOR) ---
    [p, n, e] = fileparts(atlas_path);
    temp_atlas = fullfile(out_dir, ['r' n e]); % e.g. Results/rBatlas.nii
    
    % Setup SPM Reslice
    flags.mean = false;
    flags.which = 1;      % Reslice 2nd image to match 1st
    flags.interp = 0;     % <--- 0 = NEAREST NEIGHBOR. CRITICAL!
    flags.prefix = 'r';
    
    spm_reslice({jac_path, atlas_path}, flags);
    
    % Move the result
    source_reslice = fullfile(fileparts(atlas_path), ['r' n e]);
    if exist(source_reslice, 'file')
        movefile(source_reslice, temp_atlas);
    elseif exist(fullfile(pwd, ['r' n e]), 'file')
        movefile(fullfile(pwd, ['r' n e]), temp_atlas);
    end

    % --- STEP 2: LOAD DATA ---
    V_jac = spm_vol(jac_path);
    Y_jac = spm_read_vols(V_jac);
    
    V_atl = spm_vol(temp_atlas);
    Y_atl = spm_read_vols(V_atl);
    
    % --- STEP 3: SANITIZE REGIONS ---
    % Force rounding to kill any tiny floating point errors
    Y_atl = round(Y_atl);
    
    region_ids = unique(Y_atl(:));
    region_ids(region_ids <= 0 | isnan(region_ids)) = []; 
    
    fprintf('   -> Found %d unique regions. (Should be 24)\n', length(region_ids));
    
    % --- STEP 4: CALCULATE ---
    mean_vals = zeros(length(region_ids), 1);
    std_vals  = zeros(length(region_ids), 1);
    percent_change = zeros(length(region_ids), 1);
    
    for i = 1:length(region_ids)
        id = region_ids(i);
        mask = (Y_atl == id);
        
        vals = Y_jac(mask);
        vals(isnan(vals)) = [];
        
        if ~isempty(vals)
            mean_vals(i) = mean(vals);
            std_vals(i)  = std(vals);
            % Convert 1.05 to +5.0%
            percent_change(i) = (mean_vals(i) - 1) * 100;
        end
    end
    
    % --- STEP 5: SAVE ---
    T = table(region_ids, mean_vals, std_vals, percent_change, ...
        'VariableNames', {'Region_ID', 'Mean_Jacobian', 'Std_Dev', 'Percent_Change'});
    
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end
    out_name = fullfile(out_dir, ['Stats_' bat_name '.csv']);
    writetable(T, out_name);
    
    if exist(temp_atlas, 'file'), delete(temp_atlas); end
    
    fprintf('   -> Done! Stats saved to: %s\n', out_name);
    fprintf('---------------------------------------------------\n');
end