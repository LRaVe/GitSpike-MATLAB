%% Script to generate a dataset in labeled line population 
% Date: June 2026 
% Author : Laure WOLFF

function CellMatrix = generate_and_plot_raster_ll(num_stimuli, num_repetitions, ...
    num_indi, num_neurons, t1, t2, base_rate, refrac, showing, plotting, other_figs)
    
    % Input parameter consistency check
    if num_indi > num_neurons
        error('Error: The number of the individual coding neurons cannot exceed the total num_neurons!');
    end
    
    num_trials = num_stimuli * num_repetitions;
    CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);
    
    %% 1. Generation of the Labeled Line Preference Matrix (GENERIC & AUTOMATIC)
    if num_indi == 4 && num_stimuli == 4 
        % Configuration historique de référence
        pref_matrix = [
            1, 0, 1, 0;
            1, 1, 0, 0;  % Neuron 1: Sensitive to S1, S2
            0, 0, 1, 1;  % Neuron 2: Sensitive to S3, S4
              % Neuron 3: Sensitive to S1, S3
            0, 1, 0, 1   % Neuron 4: Sensitive to S2, S4
        ];
    else
        % Génération géometrique automatique pour TOUTES les autres tailles (5x5, 8x10, etc.)
        % Assure que chaque neurone individuel possède une sélectivité propre et structurée.
        pref_matrix = zeros(num_indi, num_stimuli);
        for c_idx = 1:num_indi
            % On distribue les motifs d'activation de manière glissante/circulaire
            shift_idx = mod(c_idx - 1, num_stimuli) + 1;
            pref_matrix(c_idx, shift_idx) = 1;
            % Optionnel : donner une deuxième sensibilité si assez de stimuli pour croiser les fonctionnalités
            if num_stimuli > 2
                shift_idx_2 = mod(c_idx, num_stimuli) + 1;
                pref_matrix(c_idx, shift_idx_2) = 1;
            end
        end
    end
    
    if showing == true
        fprintf('\n--- Labeled Line Preference Matrix (Neurons %dx%d Stimuli) ---\n', num_indi, num_stimuli);
        disp(pref_matrix);
        fprintf('----------------------------------------------------------\n\n');
    end
    
    %% 2. Spikes Generation for CODING Neurons
    for c_idx = 1:num_indi
        for st = 1:num_stimuli
            
            % Attribution de taux distincts pour que la métrique de distance
            % puisse séparer les stimuli, même au sein de la préférence.
            if pref_matrix(c_idx, st) == 1
                % Stimulus préféré principal vs secondaire
                if mod(st, 2) == 1
                    local_rate = base_rate * 1.5; % Signal fort type A (ex: 45 Hz)
                else
                    local_rate = base_rate * 1.0; % Signal fort type B (ex: 30 Hz)
                end
            else
                % Bruit de fond structurel mais stable (évite le silence aléatoire)
                local_rate = base_rate * 0.1; % Baseline faible et régulière
            end

            % %% Pour tester un dataset à très fort contrast
            % if pref_matrix(c_idx, st) == 1
            %     local_rate = base_rate; % Activité intense et stable (ex: 24 Hz)
            % else
            %     local_rate = 0.001;           % Quasiment AUCUN spike (Silence radio)
            % end

            % if pref_matrix(c_idx, st) == 1
            %     local_rate = 30;  % Signal fort (30 spikes/s en moyenne)
            % else
            %     local_rate = 1;  % Bruit de fond (10 spikes/s en moyenne)
            % end
            
            for rp = 1:num_repetitions
                approx_spikes = round((t2 - t1) * local_rate * 3) + 10;
                uniform_samples = rand(1, approx_spikes);
                intervals = refrac - log(1 - uniform_samples) / local_rate;
                spikes = cumsum(intervals);
                CellMatrix{c_idx, st, rp} = spikes(spikes >= t1 & spikes <= t2);
            end
        end
    end
    
    %% 3. Generation of NON-CODING background neurons
    for st = 1:num_stimuli
        for rp = 1:num_repetitions
            for c_idx = (num_indi + 1):num_neurons
                approx_spikes = round((t2 - t1) * base_rate * 1.0) + 10;
                uniform_samples = rand(1, approx_spikes);
                intervals = refrac - log(1 - uniform_samples) / base_rate;
                spikes_noise = cumsum(intervals);
                CellMatrix{c_idx, st, rp} = spikes_noise(spikes_noise >= t1 & spikes_noise <= t2);
            end
        end
    end
    
    %% =========================================================================
    %% PLOTTING SECTION (Layout completely dynamic)
    %% =========================================================================
    if plotting == true
        if other_figs == true
            figure('Name', 'Individual Neuronal Raster Plots (LL Mode)', 'Color', 'w', 'Position', [100, 100, 1400, 900]);
            
            cols = ceil(sqrt(num_neurons * 1.25)); 
            rows = ceil(num_neurons / cols);
            
            %% Loop over each individual neuron
            for n = 1:num_neurons
                subplot(rows, cols, n);
                hold on;
                
                % Background blocks colors separation
                for st = 1:num_stimuli
                    y_start = (st - 1) * num_repetitions + 0.5;
                    y_end = st * num_repetitions + 0.5;
                    
                    if mod(st, 2) == 1
                        fill([t1 t2 t2 t1], [y_start y_start y_end y_end], [0.96 0.96 0.96], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                    end
                end
                
                % Plotting spikes lines
                t_counter = 1;
                for st = 1:num_stimuli
                    for rp = 1:num_repetitions
                        spikes = CellMatrix{n, st, rp};
                        if ~isempty(spikes)
                            for sp = 1:length(spikes)
                                line([spikes(sp), spikes(sp)], [t_counter - 0.4, t_counter + 0.4], ...
                                     'Color', 'k', 'LineWidth', 1.0);
                            end
                        end
                        t_counter = t_counter + 1;
                    end
                end
                
                % Horizontal stimulus delimitation lines
                for st_sep = 1:(num_stimuli-1)
                    sep_line = st_sep * num_repetitions + 0.5;
                    line([t1, t2], [sep_line, sep_line], 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'LineStyle', '--');
                end
                
                box on; 
                xlim([t1, t2]); 
                ylim([0.5, num_trials + 0.5]);
                set(gca, 'YDir', 'reverse');
                
                % Dynamic Tick labels handling
                if num_stimuli <= 8 % Augmenté pour accepter plus de stimuli de manière lisible
                    y_ticks = (0.5 : num_repetitions : num_trials) + (num_repetitions / 2);
                    y_labels = cell(1, num_stimuli);
                    for s = 1:num_stimuli
                        y_labels{s} = sprintf('S%d', s);
                    end
                    set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels, 'FontSize', 7);
                else
                    set(gca, 'YTick', [1, num_trials], 'YTickLabel', {'1', num2str(num_trials)}, 'FontSize', 7);
                end
                
                if n <= num_indi
                    title(sprintf('Neuron %d (Coding)', n), 'Color', [0.85 0 0], 'FontSize', 9, 'FontWeight', 'bold');
                else
                    title(sprintf('Neuron %d (Noise)', n), 'Color', [0 0.45 0.74], 'FontSize', 9, 'FontWeight', 'normal');
                end
                
                if mod(n-1, cols) == 0
                    ylabel('Stimuli / Trials', 'FontSize', 7);
                end
                if n > (num_neurons - cols)
                    xlabel('Time', 'FontSize', 7);
                end
            end
            hold off;
        end
        
        %% Multi-Subplot layout per Trial (Dynamic Fig 1 of 2018 Paper)
        figure('Name', 'Raster plot Fig 1 of the 2018s paper', 'Color', 'w');
        idx_subplot = 0; 
        for st_select = 1:num_stimuli
            for rp_select = 1:num_repetitions
                idx_subplot = idx_subplot + 1;
                subplot(num_trials, 1, idx_subplot);
                hold on;
                
                max_indi_size = 0; max_noise_size = 0;
                for c_idx = 1:num_neurons
                    sz = length(CellMatrix{c_idx, st_select, rp_select});
                    if c_idx <= num_indi
                        max_indi_size = max_indi_size + sz;
                    else
                        max_noise_size = max_noise_size + sz; 
                    end
                end
                all_indi_spikes = zeros(1, max_indi_size);
                all_noise_spikes = zeros(1, max_noise_size);
                
                ptr_i = 1; ptr_nc = 1;
                for c_idx = 1:num_neurons
                    spikes = CellMatrix{c_idx, st_select, rp_select};
                    y_pos = (num_neurons - c_idx + 1) + 4;
                    
                    if c_idx <= num_indi
                        current_color = [1 0 0];
                        len = length(spikes);
                        if len > 0
                            all_indi_spikes(ptr_i : ptr_i + len - 1) = spikes;
                            ptr_i = ptr_i + len;
                        end
                    else
                        current_color = [0.0000 0.4470 0.7410];
                        len = length(spikes);
                        if len > 0
                            all_noise_spikes(ptr_nc : ptr_nc + len - 1) = spikes;
                            ptr_nc = ptr_nc + len;
                        end
                    end
                    if ~isempty(spikes)
                        for sp = 1:length(spikes)
                            line([spikes(sp), spikes(sp)], [y_pos - 0.35, y_pos + 0.35], ...
                                 'Color', current_color, 'LineWidth', 0.8);
                        end
                    end
                end
                
                all_indi_spikes  = unique(sort(all_indi_spikes(1:ptr_i-1)));
                all_noise_spikes = unique(sort(all_noise_spikes(1:ptr_nc-1)));
                all_total_spikes = unique(sort([all_indi_spikes, all_noise_spikes]));
                
                for sp = 1:length(all_total_spikes)
                    line([all_total_spikes(sp), all_total_spikes(sp)], [1 - 0.4, 1 + 0.4], 'Color', [0 0 0], 'LineWidth', 1.3);
                end
                for sp = 1:length(all_noise_spikes)
                    line([all_noise_spikes(sp), all_noise_spikes(sp)], [2 - 0.3, 2 + 0.3], 'Color', [0.0000 0.4470 0.7410], 'LineWidth', 1.0);
                end
                for sp = 1:length(all_indi_spikes)
                    line([all_indi_spikes(sp), all_indi_spikes(sp)], [3 - 0.3, 3 + 0.3], 'Color', [1 0 0], 'LineWidth', 1.5);
                end
               
                line([t1, t2], [4.5, 4.5], 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
                line([t1, t2], [num_neurons + 4.5 - num_indi, num_neurons + 4.5 - num_indi], 'Color', [0.7 0.7 0.7], 'LineStyle', ':');
                
                box on; xlim([t1, t2]); ylim([0.5, num_neurons + 5.5]);
                
                y_ticks = [1, 2, 3, 4, num_neurons+3-num_indi, num_neurons+3];
                y_labels = {'Total', '\Sigma NC', '\Sigma Indi', num2str(num_neurons), num2str(num_indi+1), '1'};
                
                [y_ticks, unique_idx] = unique(y_ticks);
                y_labels = y_labels(unique_idx);
                
                set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels, 'FontSize', 7);
                title(sprintf('Trial : S%d-R%d', st_select, rp_select), 'FontSize', 8, 'FontWeight', 'bold');
                set(gca, 'XTickLabel', []);
            end
        end
        xlabel('Time (au)', 'FontSize', 9, 'FontWeight', 'bold');
    end
end