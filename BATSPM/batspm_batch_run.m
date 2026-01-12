function batspm_batch_run(input_folder, output_root, run_stats)
% BATSPM_BATCH_RUN
% Automates the pipeline for an entire folder.
%
% USAGE:
%   batspm_batch_run('RawData', 'Results')        -> Default (No Stats)
%   batspm_batch_run('RawData', 'Results', true)  -> With Atlas Stats

    if nargin < 3, run_stats = false; end

    % 1. Setup
    if ~exist(input_folder, 'dir'), error('Input folder missing'); end
    if ~exist(output_root, 'dir'), mkdir(output_root); end

    file_list = dir(fullfile(input_folder, '*.nii'));
    num_files = length(file_list);
    
    fprintf('---------------------------------------------------\n');
    fprintf('BATSPM BATCH MODE (%d subjects)\n', num_files);
    fprintf('Atlas Statistics: %s\n', string(run_stats));
    fprintf('---------------------------------------------------\n');
    
    % 2. Loop
    for i = 1:num_files
        filename = file_list(i).name;
        full_path = fullfile(input_folder, filename);
        
        fprintf('\n>>> Processing %d/%d: %s\n', i, num_files, filename);
        try
            tic;
            % Pass the 'run_stats' flag down to the main runner
            batspm_run(full_path, output_root, run_stats);
            fprintf('   [OK] %.1fs\n', toc);
        catch ME
            fprintf('   [ERROR] %s failed: %s\n', filename, ME.message);
        end
    end
    fprintf('\nBatch Complete.\n');
end