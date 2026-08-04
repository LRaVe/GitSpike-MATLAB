%% Script to use the user's data
% Date: June 2026
% Author : Laure WOLFF

function spikes = import_spikes_from_txt(filename, num_stimuli, num_repetitions, num_neurons)
% IMPORT_SPIKES_FROM_TXT Dynamic data loading from text file to a 3D Cell Matrix.
%
%   spikes = IMPORT_SPIKES_FROM_TXT(filename, num_stimuli, num_repetitions, num_neurons) 
%   reads a structured text file line by line and reconstructs the multi-neuron 
%   spike train dataset.
%
%   INPUTS:
%       data/file parameters:
%       - filename        : Character vector or string specifying the path to the text file.
%       - num_stimuli     : Exact number of stimuli (S).
%       - num_repetitions : Exact number of repetitions/trials per stimulus (R).
%       - num_neurons     : Exact number of recorded neurons (N).
%
%   OUTPUTS:
%       - spikes          : 3D cell array of size {num_neurons x num_stimuli x num_repetitions} 
%                           formatted and ready for downstream analysis algorithms.
%
%   COMPATIBILITY:
%       - Perfectly mirrors the writing format produced by export_spikes_to_txt.m.
%
%   Author: Laure WOLFF
%   Date: June 2026

    % 1. Read all lines from the file
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
            parsed_nums = str2num(cleaned_line); %#ok<ST2NM> Deliberatly using rather str2double
            
            if isempty(parsed_nums)
                raw_lines{total_lines} = []; 
            else
                raw_lines{total_lines} = parsed_nums(:)'; % Force horizontal row vector
            end 
        end
    end
    fclose(fid);
    raw_lines = raw_lines(1:total_lines); % Adjust to true size
       
    % 3. Protocol Validation Check (Sécurité dimensions)
    expected_lines = num_stimuli * num_repetitions * num_neurons;  % {S,R,N}
    if total_lines ~= expected_lines
        error(['Fatal dimension mismatch error!\n' ...
               'According to your inputs (%d Stimuli x %d Reps x %d Neurons), the file should have %d lines.\n' ...
               'However, the file contains %d lines.'], ...
               num_stimuli, num_repetitions, num_neurons, expected_lines, total_lines);
    end
        
    % 4. Structured distribution into a temporary 2D CellMatrix {Stimulus, Repetition}
    temp = cell(num_stimuli, num_repetitions);
    line_pointer = 1;
    
    for s = 1:num_stimuli
        for r = 1:num_repetitions
            trial_neurons_block = cell(num_neurons, 1);
            for n = 1:num_neurons
                trial_neurons_block{n} = raw_lines{line_pointer};
                line_pointer = line_pointer + 1;
            end
            temp{s, r} = trial_neurons_block;
        end
    end
    
    % 5. CONVERSION BLOCK: Re-shaping into the definitive 3D matrix {N, S, R}
    spikes = cell(num_neurons, num_stimuli, num_repetitions);
    for s = 1:num_stimuli
        for r = 1:num_repetitions
            for n = 1:num_neurons
                spikes{n, s, r} = temp{s, r}{n};
            end
        end
    end
    
    fprintf('spikes {%dx%dx%d} successfully generated for algorithms.\n', num_neurons, num_stimuli, num_repetitions);
end


