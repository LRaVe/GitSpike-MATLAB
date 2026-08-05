function export_spikes_to_txt(spikes,filename)
% EXPORT_SPIKES_TO_TXT Export a spike dataset to a formatted text file.
%
%   EXPORT_SPIKES_TO_TXT(spikes, filename) takes a 3D cell array 
%   {Neuron, Stimulus, Repetition} containing spike timestamps and exports 
%   them into a structured text file.
%
%   INPUTS:
%       - spikes   : 3D cell array of dimensions [num_neurons x num_stimuli x num_repetitions] 
%                    holding the spike timestamp vectors for each trial.
%       - filename : Character vector or string specifying the destination file path.
%
%   FILE FORMAT:
%       - Each line in the generated text file corresponds to the spike train 
%         of a single neuron for a single trial.
%       - Timestamps are formatted with high precision ('%.17e') and separated by spaces.
%       - Loop nesting order: Stimulus -> Repetition -> Neuron.
%
%   COMPATIBILITY:
%       - Fully compatible with the complementary function import_spikes_from_txt.m.

    fid = fopen(filename,'w');

    if fid == -1
        error('Could not create file "%s".',filename);
    end

    [N,S,R] = size(spikes);

    fprintf('=== EXPORTING SPIKE DATASET ===\n');

    for s = 1:S

        for r = 1:R

            for n = 1:N

                train = spikes{n,s,r};

                if isempty(train)

                    fprintf(fid,'\n');

                else

                    train = sort(train(:));

                    fprintf(fid,'%.17e ',train(1:end-1));
                    fprintf(fid,'%.17e\n',train(end));

                end

            end

        end

    end

    fclose(fid);

    fprintf('File "%s" successfully created.\n',filename);
    fprintf('Exported %d neurons x %d stimuli x %d repetitions (%d lines).\n',...
        N,S,R,N*S*R);

end