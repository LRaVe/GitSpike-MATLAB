%% Script to generate a dataset in labeled line population 
% Date: June 2026 
% Author : Laure WOLFF (Optimized Raster Plot Vectorization)
function CellMatrix = generate_and_plot_raster_ll(num_stimuli, num_repetitions, ...
    num_indi, num_neurons, t1, t2, base_rate, refrac, jitter_std, showing, plotting, other_figs)

    % Input parameter consistency check
    if num_indi > num_neurons
        error('Error: The number of individual coding neurons cannot exceed total num_neurons!');
    end
    
    num_trials = num_stimuli * num_repetitions;
    CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);
    
    %% 1. Generation of the Labeled Line Preference Matrix
    if num_indi == 4 && num_stimuli == 4 
        pref_matrix = [
            1, 1, 0, 0;  % Neuron 1: Sensitive to S1, S2
            0, 0, 1, 1;  % Neuron 2: Sensitive to S3, S4
            1, 0, 1, 0;  % Neuron 3: Sensitive to S1, S3
            0, 1, 0, 1   % Neuron 4: Sensitive to S2, S4
        ];
    else
        % Pure random Bernoulli distribution (No pre-determined active counts)
        pref_matrix = rand(num_indi, num_stimuli) < 0.5;
        
        % Safeguard: Ensure no silent neurons (empty rows)
        for c_idx = 1:num_indi
            if sum(pref_matrix(c_idx, :)) == 0
                pref_matrix(c_idx, randi(num_stimuli)) = 1;
            end
        end
        % Safeguard: Ensure all stimuli are encoded (empty columns)
        for s_idx = 1:num_stimuli
            if sum(pref_matrix(:, s_idx)) == 0
                pref_matrix(randi(num_indi), s_idx) = 1;
            end
        end
    end
    
    if showing == true
        fprintf('\n--- Labeled Line Preference Matrix (Neurons %dx%d Stimuli) ---\n', num_indi, num_stimuli);
        disp(double(pref_matrix));
        fprintf('----------------------------------------------------------\n\n');
    end
    
    %% 2. Spikes Generation for CODING Neurons
    for c_idx = 1:num_indi
        
        % --- STEP A: Generate shared anchor temporal signature (Baseline) ---
        approx_spikes = round((t2 - t1) * base_rate * 2) + 10;
        uniform_samples = rand(1, approx_spikes);
        intervals = refrac - log(1 - uniform_samples) / base_rate;
        baseline_spikes = t1 + cumsum(intervals);
        baseline_spikes = baseline_spikes(baseline_spikes >= t1 & baseline_spikes <= t2);
        num_spikes = length(baseline_spikes);
        
        if num_spikes == 0, num_spikes = 1; baseline_spikes = (t1 + t2) / 2; end
        
        % --- STEP B: Distribute jittered trials or low-frequency background ---
        for st = 1:num_stimuli
            if pref_matrix(c_idx, st) == 1
                % Jittered common template execution
                for rp = 1:num_repetitions
                    shifts = randn(1, num_spikes) * jitter_std;
                    jittered_train = baseline_spikes + shifts;
                    
                    jittered_train(jittered_train < t1) = t1;
                    jittered_train(jittered_train > t2) = t2;
                    CellMatrix{c_idx, st, rp} = sort(jittered_train);
                end
            else
                % Asynchronous random Poisson noise
                low_noise_rate = base_rate * 2.5; 
                approx_spikes_noise = round((t2 - t1) * low_noise_rate * 2) + 10;
                for rp = 1:num_repetitions
                    uniform_samples_noise = rand(1, approx_spikes_noise);
                    intervals_noise = refrac - log(1 - uniform_samples_noise) / low_noise_rate;
                    spikes_noise = t1 + cumsum(intervals_noise);
                    CellMatrix{c_idx, st, rp} = spikes_noise(spikes_noise >= t1 & spikes_noise <= t2);
                end
            end
        end
    end
    
    %% 3. Generation of NON-CODING background neurons
    for c_idx = (num_indi + 1):num_neurons
        for st = 1:num_stimuli
            low_noise_rate = base_rate * 2.5; 
            approx_spikes_noise = round((t2 - t1) * low_noise_rate * 2) + 10;
            for rp = 1:num_repetitions
                uniform_samples_noise = rand(1, approx_spikes_noise);
                intervals_noise = refrac - log(1 - uniform_samples_noise) / low_noise_rate;
                spikes_noise = t1 + cumsum(intervals_noise);
                CellMatrix{c_idx, st, rp} = spikes_noise(spikes_noise >= t1 & spikes_noise <= t2);
            end
        end
    end
    
    %% =========================================================================
    %% PLOTTING SECTION
    %% =========================================================================
    if plotting == true
        %% FIGURE 1: Individual Neuronal Raster Plots (LL Mode)
        if other_figs == true
            figure('Name', 'Individual Neuronal Raster Plots (LL Mode)', 'Color', 'w');
            cols = ceil(sqrt(num_neurons * 1.25)); 
            rows = ceil(num_neurons / cols);
            
            for n = 1:num_neurons
                subplot(rows, cols, n);
                hold on;
                
                % Background stimulus blocks shading
                for st = 1:num_stimuli
                    y_start = (st - 1) * num_repetitions + 0.5;
                    y_end = st * num_repetitions + 0.5;
                    if mod(st, 2) == 1
                        fill([t1 t2 t2 t1], [y_start y_start y_end y_end], [0.96 0.96 0.96], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
                    end
                end
                
                % Step 1: Pre-calculate exact spike count to avoid dynamic growing vectors
                max_spikes = 0;
                for st = 1:num_stimuli
                    for rp = 1:num_repetitions
                        max_spikes = max_spikes + length(CellMatrix{n, st, rp});
                    end
                end
                
                % Step 2: Preallocate line array segments [X1; X2; NaN]
                total_elements = max_spikes * 3;
                X_lines = NaN(1, total_elements);
                Y_lines = NaN(1, total_elements);
                idx_pointer = 1;
                
                t_counter = 1;
                % Step 3: Fast buffer filling
                for st = 1:num_stimuli
                    for rp = 1:num_repetitions
                        spikes = CellMatrix{n, st, rp};
                        num_spikes_current = length(spikes);
                        
                        if num_spikes_current > 0
                            idx_range = idx_pointer : (idx_pointer + (num_spikes_current * 3) - 1);
                            
                            X_spikes = [spikes; spikes; NaN(1, num_spikes_current)];
                            Y_spikes = [ones(1, num_spikes_current) * (t_counter - 0.4); ...
                                        ones(1, num_spikes_current) * (t_counter + 0.4); ...
                                        NaN(1, num_spikes_current)];
                                    
                            X_lines(idx_range) = X_spikes(:)';
                            Y_lines(idx_range) = Y_spikes(:)';
                            
                            idx_pointer = idx_pointer + (num_spikes_current * 3);
                        end
                        t_counter = t_counter + 1;
                    end
                end
                
                % Crop and display vectorised graphics block
                if idx_pointer > 1
                    X_lines = X_lines(1:idx_pointer-1);
                    Y_lines = Y_lines(1:idx_pointer-1);
                    line(X_lines, Y_lines, 'Color', 'k', 'LineWidth', 1.0);
                end
                
                % Dashed boundaries lines between stimuli
                for st_sep = 1:(num_stimuli-1)
                    sep_line = st_sep * num_repetitions + 0.5;
                    line([t1, t2], [sep_line, sep_line], 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'LineStyle', '--');
                end
                
                box on; xlim([t1, t2]); ylim([0.5, num_trials + 0.5]);
                set(gca, 'YDir', 'reverse');
                
                if num_stimuli <= 8
                    y_ticks = (0.5 : num_repetitions : num_trials) + (num_repetitions / 2);
                    y_labels = arrayfun(@(s) sprintf('S%d', s), 1:num_stimuli, 'UniformOutput', false);
                    set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels, 'FontSize', 7);
                else
                    set(gca, 'YTick', [1, num_trials], 'YTickLabel', {'1', num2str(num_trials)}, 'FontSize', 7);
                end
                
                if n <= num_indi
                    title(sprintf('Neuron %d (Coding)', n), 'Color', [0.85 0 0], 'FontSize', 9, 'FontWeight', 'bold');
                else
                    title(sprintf('Neuron %d (Noise)', n), 'Color', [0 0.45 0.74], 'FontSize', 9, 'FontWeight', 'normal');
                end
                if mod(n-1, cols) == 0, ylabel('Stimuli / Trials', 'FontSize', 7); end
                if n > (num_neurons - cols), xlabel('Time', 'FontSize', 7); end
            end
        end
        % %% FIGURE 2 : Échantillons de Trials (Inspiré Fig 1, Paper 2018)

        % figure('Name', 'Raster plot Fig 1 (Selected Trials Samples)', 'Color', 'w');
        % max_subplots = min(4, num_trials);
        % trials_to_plot = round(linspace(1, num_trials, max_subplots));
        %
        % for idx_sub = 1:length(trials_to_plot)
        % subplot(length(trials_to_plot), 1, idx_sub);
        % hold on;
        %
        % current_trial = trials_to_plot(idx_sub);
        % st_select = ceil(current_trial / num_repetitions);
        % rp_select = mod(current_trial - 1, num_repetitions) + 1;
        %
        % for c_idx = 1:num_neurons
        % spikes = CellMatrix{c_idx, st_select, rp_select};
        % y_pos = (num_neurons - c_idx + 1) + 4;
        %
        % current_color = [0.0000 0.4470 0.7410];
        % if c_idx <= num_indi, current_color = [1 0 0]; end
        %
        % if ~isempty(spikes)
        % X_sp = [spikes; spikes; NaN(1, length(spikes))];
        % Y_sp = [ones(1, length(spikes))*(y_pos - 0.35); ones(1, length(spikes))*(y_pos + 0.35); NaN(1, length(spikes))];
        % line(X_sp(:)', Y_sp(:)', 'Color', current_color, 'LineWidth', 0.8);
        % end
        % end
        %
        % % Sommes de populations cumulées
        % all_indi = cell2mat(reshape(CellMatrix(1:num_indi, st_select, rp_select), 1, []));
        % all_noise = cell2mat(reshape(CellMatrix(num_indi+1:end, st_select, rp_select), 1, []));
        % all_total = [all_indi, all_noise];
        %
        % if ~isempty(all_total)
            % X_t = [all_total; all_total; NaN(1, length(all_total))];
            % Y_t = [ones(1, length(all_total))*0.6; ones(1, length(all_total))*1.4; NaN(1, length(all_total))];
            % line(X_t(:)', Y_t(:)', 'Color', [0 0 0], 'LineWidth', 1.2);
        % end
        % if ~isempty(all_noise)
        % X_n = [all_noise; all_noise; NaN(1, length(all_noise))];
        % Y_n = [ones(1, length(all_noise))*1.7; ones(1, length(all_noise))*2.3; NaN(1, length(all_noise))];
        % line(X_n(:)', Y_n(:)', 'Color', [0.0000 0.4470 0.7410], 'LineWidth', 1.0);
        % end
        % if ~isempty(all_indi)
        % X_i = [all_indi; all_indi; NaN(1, length(all_indi))];
        % Y_i = [ones(1, length(all_indi))*2.7; ones(1, length(all_indi))*3.3; NaN(1, length(all_indi))];
        % line(X_i(:)', Y_i(:)', 'Color', [1 0 0], 'LineWidth', 1.2);
        % end
        %
        % line([t1, t2], [4.5, 4.5], 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
        % line([t1, t2], [num_neurons + 4.5 - num_indi, num_neurons + 4.5 - num_indi], 'Color', [0.7 0.7 0.7], 'LineStyle', ':');
        %
        % box on; xlim([t1, t2]); ylim([0.5, num_neurons + 5.5]);
        %
        % y_ticks = unique([1, 2, 3, 5, num_neurons + 4]);
        % y_labels = {'Total', '\Sigma NC', '\Sigma Indi', 'Neurone N', 'Neurone 1'};
        % set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels, 'FontSize', 7);
        %
        % title(sprintf('Trial: S%d - Rep%d', st_select, rp_select), 'FontSize', 8, 'FontWeight', 'bold');
        % if idx_sub < length(trials_to_plot), set(gca, 'XTickLabel', []); end
        % end
        % xlabel('Time (au)', 'FontSize', 9, 'FontWeight', 'bold');
    end
end