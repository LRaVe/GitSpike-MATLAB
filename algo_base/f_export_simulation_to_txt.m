function f_export_simulation_to_txt(CellMatrix, filename)
% f_export_simulation_to_txt - Version compatible 3D {N, S, R} pour ton script
% Chaque ligne du fichier correspondra STRICTEMENT au train de spikes d'un neurone pour un essai donné.

    fid = fopen(filename, 'w');
    if fid == -1
        error('Could not create or open file: %s', filename);
    end
    
    % On récupère les dimensions depuis ta matrice 3D d'origine
    [num_neurons, num_stimuli, num_repetitions] = size(CellMatrix);
    
    fprintf('=== EXPORTING SIMULATION DATA (COMPATIBLE FORMAT) ===\n');
    
    % On parcourt dans l'ordre exact : Stimuli -> Repetitions -> Neurones
    for s = 1:num_stimuli
        for r = 1:num_repetitions
            for n = 1:num_neurons
                
                spikes = CellMatrix{n, s, r};

                if isempty(spikes)
                    spikes_sorted = [];
                else
                    spikes_col = spikes(:); 
                    spikes_sorted = sort(spikes_col);
                end
                
                if ~isempty(spikes_sorted)
                    for i = 1:length(spikes_sorted)-1
                        fprintf(fid, '%.17e ', spikes_sorted(i));
                    end
                    fprintf(fid, '%.17e \n', spikes_sorted(end));
                else
                    fprintf(fid, '\n'); 
                end
                
            end
        end
    end
    
    fclose(fid);
    total_lines = num_stimuli * num_repetitions * num_neurons;
    fprintf('File "%s" successfully created! Wrote %d lines.\n', filename, total_lines);
end