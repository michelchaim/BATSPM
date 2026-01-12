function batspm_run(T2_path, out_dir, do_stats)
% BATSPM_RUN - Main pipeline wrapper for Bat MRI Analysis
% =========================================================
% USAGE:
%   batspm_run(image, output)           -> Standard run
%   batspm_run(image, output, true)     -> Run + Atlas Statistics

    % --- 0. Setup Defaults ---
    if nargin < 3, do_stats = false; end 

    % 1. Get Bat Name & Check Extension
    [~, fname, fext] = fileparts(char(T2_path));
    
    % GUARD: Reject .gz or .zip
    if strcmpi(fext, '.gz') || strcmpi(fext, '.zip')
        error('STOP! Please unzip this file manually. This script only accepts .nii files.');
    end
    
    % GUARD: Reject anything that isn't .nii
    if ~strcmpi(fext, '.nii')
        error('STOP! The input file must be a .nii file. You provided: %s', fext);
    end

    % 2. Setup Results Structure
    fname_clean = regexprep(fname, '[^a-zA-Z0-9]', '');
    bat_dir = fullfile(out_dir, fname_clean);
    subject_dir = bat_dir; 
    
    if ~exist(bat_dir, 'dir'), mkdir(bat_dir); end
    new_T2_path = fullfile(bat_dir, [fname_clean '.nii']);
    
    % 3. Copy the Data
    if ~exist(new_T2_path, 'file')
        copyfile(T2_path, new_T2_path);
        fprintf('Moved data to workspace: %s\n', bat_dir);
    end
    
    % --- PIPELINE START ---
    
    % Locate Toolbox & TPM
    batspm_dir = fileparts(mfilename('fullpath'));
    
    % NOTE: Checking for 'Bat_tpm.nii' as per your code. 
    % If you are using standard SPM TPM, change this to 'TPM.nii'
    tpm_path = fullfile(batspm_dir, 'Bat_tpm.nii'); 
    
    if ~exist(tpm_path, 'file')
        % Fallback: Check for standard name if Bat_tpm is missing
        tpm_path = fullfile(batspm_dir, 'TPM.nii');
        if ~exist(tpm_path, 'file')
             error('TPM not found! Checked for Bat_tpm.nii and TPM.nii in %s', batspm_dir); 
        end
    end

    % 4. Reslice to Isotropic (The "Pancake" Fix)
    if exist('batspm_reslice_iso', 'file')
        batspm_reslice_iso(new_T2_path, 0.4);
    else
        warning('batspm_reslice_iso not found! Skipping geometry fix.');
    end

    % 5. Run Segmentation
    fprintf('Running Bat Segmentation...\n');
    seg = batspm_segment(new_T2_path, tpm_path);
    
    % 6. Normalization 
    norm_res = [];
    if exist('batspm_normalize', 'file')
        norm_res = batspm_normalize(seg, subject_dir); 
    end
        
    % 7. Jacobian 
    jac = [];
    if exist('batspm_jacobian', 'file') && ~isempty(norm_res)
        jac = batspm_jacobian(seg, subject_dir);
    
        % 7b. Create QC Snapshot (Safe Version)
        if exist('batspm_save_snapshot', 'file') && ~isempty(norm_res) && ~isempty(jac)
             batspm_save_snapshot(norm_res.wT2, jac.det, subject_dir, fname_clean);
        end
        
        % 8. Atlas statistics
        if do_stats
            atlas_path = fullfile(batspm_dir, 'Batlas_labels.nii');
            if exist('batspm_atlas_stats', 'file') && exist(atlas_path, 'file')
                batspm_atlas_stats(jac.det, atlas_path, subject_dir, fname_clean);
            else
                warning('Skipping Stats: Missing batspm_atlas_stats.m or Batlas.nii');
            end
        end
    end
    
    % 9. Regional Volumes
    if exist('batspm_extract_volumes', 'file')
        volumes = batspm_extract_volumes(seg, subject_dir);
        
        % 10. Save Summary 
        if exist('batspm_save_results', 'file')
            batspm_save_results(volumes, fname_clean, subject_dir);
        end
    end
    
    fprintf('Pipeline complete for %s.\n', fname_clean);
end
% =========================================================
% Sub-function: Segmentation
% =========================================================
function seg = batspm_segment(T2_path, tpm_path)
    
    spm('Defaults','fMRI');
    spm_jobman('initcfg');

    % --- SAFETY CHECK ---
    fprintf('DEBUG: Looking for file at: %s\n', T2_path);
    if exist(T2_path, 'file') ~= 2
        error('STOP! The file does not exist on the disk here: %s', T2_path);
    else
        fprintf('DEBUG: File found! Proceeding to SPM...\n');
    end
    
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.preproc.channel.vols    = {T2_path}; % Clean path
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.001; 
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60; 
    matlabbatch{1}.spm.spatial.preproc.channel.write   = [0 1]; 
    
    % --- Tissue Setup ---
    for k = 1:3
        matlabbatch{1}.spm.spatial.preproc.tissue(k).tpm    = {[tpm_path ',' num2str(k)]};
        matlabbatch{1}.spm.spatial.preproc.tissue(k).ngaus  = 2;
        matlabbatch{1}.spm.spatial.preproc.tissue(k).native = [1 0]; 
        matlabbatch{1}.spm.spatial.preproc.tissue(k).warped = [0 0]; 
    end
    
    % --- WARPING ---
    matlabbatch{1}.spm.spatial.preproc.warp.mrf     = 1; 
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 1; 
    matlabbatch{1}.spm.spatial.preproc.warp.reg     = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg  = ''; 
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm    = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp    = 1; 
    matlabbatch{1}.spm.spatial.preproc.warp.write   = [1 1]; 
    
    spm_jobman('run', matlabbatch);
    
    [pth, nam, ext] = fileparts(T2_path);
    seg.T2           = T2_path;
    seg.biascorr     = fullfile(pth, ['m' nam ext]);
    seg.gm           = fullfile(pth, ['c1' nam ext]);
    seg.wm           = fullfile(pth, ['c2' nam ext]);
    seg.csf          = fullfile(pth, ['c3' nam ext]);
    seg.defforward   = fullfile(pth, ['y_' nam ext]);
    seg.definverse   = fullfile(pth, ['iy_' nam ext]);
end

% =========================================================
% Sub-function: Normalization
% =========================================================
function norm_res = batspm_normalize(seg, out_dir)
    fprintf('Running Normalization...\n');
    spm('Defaults','fMRI');
    spm_jobman('initcfg');
    
    def_field = seg.defforward; 
    
    % Warp T2, GM, WM, and CSF
    images_to_warp = {
        seg.T2;
        seg.gm;
        seg.wm;
        seg.csf; 
    };
    matlabbatch = [];
    matlabbatch{1}.spm.spatial.normalise.write.subj.def      = {def_field};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample = images_to_warp;
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox  = [0.5 0.5 0.5]; % Bat resolution
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb   = NaN(2,3);
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 4;
    
    spm_jobman('run', matlabbatch);
    
    [p, n, e] = fileparts(seg.T2);
    norm_res.wT2 = fullfile(p, ['w' n e]);
    [p, n, e] = fileparts(seg.gm);
    norm_res.wgm = fullfile(p, ['w' n e]);
    [p, n, e] = fileparts(seg.wm);
    norm_res.wwm = fullfile(p, ['w' n e]);
    
    fprintf('Normalization complete.\n');
end

% =========================================================
% Sub-function: Jacobian (Final Force-Fix Version)
% =========================================================
function jac = batspm_jacobian(seg, out_dir) 
    fprintf('Calculating Jacobian...\n');
    
    def_field = seg.defforward; 
    [def_path, def_name, def_ext] = fileparts(def_field);
    
    % --- BUILD JOB ---
    job = struct();
    job.comp{1}.def = {def_field};
    
    % 1. Pass the ORIGINAL name (without 'j_') 
    % SPM will verify the input and add the 'j_' prefix automatically.
    % This fixes the "j_j_" issue.
    job.out{1}.savejac.ofname = [def_name def_ext]; 
    job.out{1}.savejac.preserve = 0; 
    job.out{1}.savejac.savedir = {def_path}; % We still try to set this
    
    % --- EXECUTE ---
    try
        if exist('spm_deformations', 'file')
             spm_deformations(job); % SPM25
        elseif exist('spm_defs', 'file')
             % Translate for SPM12
             job.out{1} = rmfield(job.out{1}, 'savejac');
             job.out{1}.jac.preserve = 0;
             spm_defs(job);
        end
    catch ME
        warning('SPM calculation warning: %s', ME.message);
    end
    
    % --- THE CLEANUP (Force Move & Rename) ---
    % 1. Define what we WANT (The perfect path in your Results folder)
    final_name = ['j_' def_name def_ext];
    final_path = fullfile(out_dir, final_name);
    
    % 2. Define what SPM likely CREATED (In current folder or def folder)
    % We check both possible locations where SPM might have dropped it.
    possible_locs = {
        fullfile(def_path, final_name), ...      % Inside the T2 folder
        fullfile(pwd, final_name), ...           % In the current BATSPM folder
        fullfile(pwd, ['j_j_' def_name def_ext]) % In case it double-named it again
    };
    
    found = false;
    for i = 1:length(possible_locs)
        if exist(possible_locs{i}, 'file')
            % We found it! Move it to the right place.
            movefile(possible_locs{i}, final_path);
            fprintf('Moved Jacobian to results: %s\n', final_path);
            found = true;
            break;
        end
    end
    
    if found
        jac.det = final_path;
    else
        warning('Could not locate the Jacobian file. Checked current and results folders.');
        jac = [];
    end
end

% =========================================================
% Sub-function: Save Results 
% =========================================================
function batspm_save_results(volumes, bat_name, out_dir)
    csv_file = fullfile(out_dir, [bat_name '_volumes.csv']);
    T = table({bat_name}, volumes.GM, volumes.WM, volumes.CSF, volumes.TIV, ...
        'VariableNames', {'ID', 'GM_mm3', 'WM_mm3', 'CSF_mm3', 'TIV_mm3'});
    writetable(T, csv_file);
    fprintf('Results saved to: %s\n', csv_file);
end