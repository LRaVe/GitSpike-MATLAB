function f_export_simulation_to_txt(CellMatrix, filename)
% f_export_simulation_to_txt - Exports the simulated CellMatrix to a flat text file.
% Each line in the file will correspond to a single trial, containing all sorted spikes.
%
% Author: Laure WOLFF - June 2026

    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not create or open file: %s', filename);
    end
    
    [num_stimuli, num_repetitions] = size(CellMatrix);
    
    fprintf('=== EXPORTING SIMULATION DATA ===\n');
    
    % Iterate through the CellMatrix trial by trial (translating to line by line in the file)
    for s = 1:num_stimuli
        for r = 1:num_repetitions
            
            % 1. Retrieve the data for the current trial
            trial_data = CellMatrix{s, r};
            all_spikes = [];
            
            % 2. Gather spikes from ALL neurons for this specific trial
            if iscell(trial_data)
                % If the CellMatrix contains a cell array of neurons {N, 1}
                for n = 1:length(trial_data)
                    all_spikes = [all_spikes; trial_data{n}(:)];
                end
            else
                % If it is a standard numeric array
                all_spikes = trial_data(:);
            end
            
            % 3. Crucial step: Sort all gathered spikes in chronological order
            all_spikes_sorted = sort(all_spikes);
            
            % 4. Write these values onto a single line, separated by a space
            if ~isempty(all_spikes_sorted)
                for i = 1:length(all_spikes_sorted)-1
                    fprintf(fid, '%.3f ', all_spikes_sorted(i));
                end
                % For the last element, append the newline character (\n)
                fprintf(fid, '%.3f\n', all_spikes_sorted(end));
            else
                % If the trial contains no spikes, still append a newline to keep file structure
                fprintf(fid, '\n');
            end
        end
    end
    
    fclose(fid);
    fprintf('File "%s" successfully created! Wrote %d lines.\n', filename, num_stimuli * num_repetitions);
end