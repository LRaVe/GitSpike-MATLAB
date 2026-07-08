function export_spikes_to_txt(spikes,filename)
% EXPORT_SPIKES_TO_TXT
%
% Export a spike dataset stored as a cell array {Neuron,Stimulus,Repetition}
% into a text file.
%
% Each line of the text file corresponds to one neuron's spike train
% for one trial.
%
% Order:
%   Stimulus -> Repetition -> Neuron
%
% This file is fully compatible with import_spikes_from_txt.m

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