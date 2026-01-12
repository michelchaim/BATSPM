function batspm_save_snapshot(anat_path, jac_path, out_dir, bat_name)
% BATSPM_SAVE_SNAPSHOT (Standard Version)
% - Anatomy: Grayscale (TrueColor fixed)
% - Jacobian: Jet Colormap
% - Sensitivity: Standard (0.5 to 1.5) - Shows only major changes.

    try
        fprintf('Generating QC Snapshot...\n');
        
        % 1. Load Data
        Va = spm_vol(anat_path);
        Ya = spm_read_vols(Va);
        Vj = spm_vol(jac_path);
        Yj = spm_read_vols(Vj);
        
        % 2. Pick Center Slice
        z_idx = round(size(Ya, 3) / 2);
        slice_anat = rot90(Ya(:, :, z_idx));
        slice_jac  = rot90(Yj(:, :, z_idx));
        
        % 3. Setup Figure
        f = figure('Visible', 'off', 'Color', 'k', 'Position', [100, 100, 1000, 500]);
        
        % --- LEFT PANEL: ANATOMY (Static Grayscale) ---
        ax1 = subplot(1, 2, 1);
        
        % Normalize to 0-1 and stack to RGB so colors "freeze"
        anat_norm = slice_anat - min(slice_anat(:));
        if max(anat_norm(:)) > 0
            anat_norm = anat_norm / max(anat_norm(:));
        end
        anat_rgb = repmat(anat_norm, [1 1 3]); 
        
        image(ax1, anat_rgb); 
        axis(ax1, 'image', 'off');
        title(ax1, 'Normalized Anatomy', 'Color', 'w');
        
        % --- RIGHT PANEL: JACOBIAN (Rainbow) ---
        ax2 = subplot(1, 2, 2);
        imagesc(ax2, slice_jac);
        
        % Apply Rainbow map
        colormap(ax2, 'jet');      
        colorbar(ax2, 'Color', 'w');
        
        % *** STANDARD CONTRAST (0.5 to 1.5) ***
        % This is the original, less sensitive scale.
        % 1.0 = Green. 1.5 = Red (Huge Growth). 0.5 = Blue (Huge Shrink).
        caxis(ax2, [0.5 1.5]); 
        
        axis(ax2, 'image', 'off');
        title(ax2, 'Jacobian Determinant', 'Color', 'w');
        
        sgtitle(['QC: ' bat_name], 'Color', 'w');
        
        % 4. Save
        snap_name = fullfile(out_dir, ['QC_' bat_name '.jpg']);
        saveas(f, snap_name);
        close(f);
        
        fprintf('Snapshot saved: %s\n', snap_name);
        
    catch ME
        warning('Snapshot creation failed: %s', ME.message);
        if exist('f', 'var'), close(f); end
    end
end