%% Script to use the user's data
% Date: June 2026
% Author : Laure WOLFF

function [CellMatrix, num_neurons] = f_import_data_secure(filename, num_stimuli, num_repetitions)
% f_import_data_secure - Dynamic data loading that resolves mathematical ambiguity.
%
% Inputs:
%   filename        : Text file name (e.g., 'Donnees_Simulees.txt')
%   num_stimuli     : Exact number of Stimuli (S)
%   num_repetitions : Exact number of Repetitions (R)
%
% Outputs:
%   CellMatrix      : Cell matrix of size {S, R} ready for algorithms
%   num_neurons     : Calculated and verified number of neurons

    % 1. Read all lines from the file and count the total
    fid = fopen(filename, 'r');
    if fid == -1
        error('Could not open file: %s', filename);
    end
    
    raw_lines = cell(1, 2000); % Pre-allocation
    total_lines = 0;
    
    while ~feof(fid)
        tline = fgetl(fid);
        if ischar(tline) && ~isempty(strtrim(tline))
            total_lines = total_lines + 1;
            cleaned_line = strtrim(tline); 
            parsed_nums = str2num(cleaned_line); %#ok<ST2NM> , str2num deberately used 
            if isempty(parsed_nums)
                raw_lines{total_lines} = []; 
            else
                raw_lines{total_lines} = parsed_nums(:)'; % On force une ligne horizontale
            end 
        end
    end
    fclose(fid);
    raw_lines = raw_lines(1:total_lines); % Adjust to true size

    % 2. Mathematical safety check: calculate trials needed per neuron
    trials_per_neuron = num_stimuli * num_repetitions;
    
    if rem(total_lines, trials_per_neuron) ~= 0
        error(['Fatal dimension mismatch error!\nThe file contains %d lines.\n' ...
               'With your protocol (%d Stimuli x %d Reps = %d trials), ' ...
               'the result would be %.2f neurons, which is biologically impossible.'], ...
               total_lines, num_stimuli, num_repetitions, trials_per_neuron, total_lines / trials_per_neuron);
    end
    
    % Determine the exact number of neurons
    num_neurons = total_lines / trials_per_neuron;
    
    fprintf('=== EXPERIMENTAL PROTOCOL DETECTED ===\n');
    fprintf('-> Number of stimuli (S)    : %d\n', num_stimuli);
    fprintf('-> Number of repetitions (R): %d\n', num_repetitions);
    fprintf('-> Number of neurones (N)   : %d\n\n', num_neurons);

    % 3. Structured distribution into CellMatrix {Stimulus, Repetition}
    CellMatrix = cell(num_stimuli, num_repetitions);
    
    % Use a pointer to track which line of the text file is being read
    line_pointer = 1;
    
    for s = 1:num_stimuli
        for r = 1:num_repetitions
            
            % For each cell {s, r}, store the spike trains of our N neurons
            trial_neurons_block = cell(num_neurons, 1);
            
            for n = 1:num_neurons
                % Extract the current line of the file for neuron n
                trial_neurons_block{n} = raw_lines{line_pointer};
                
                % Advance by one line in the text file
                line_pointer = line_pointer + 1;
            end
            
            % Assign this block of N neurons to the correct experimental condition
            CellMatrix{s, r} = trial_neurons_block;
        end
    end
    
    fprintf('CellMatrix {%dx%d} successfully generated for %d neurons.\n', num_stimuli, num_repetitions, num_neurons);
end