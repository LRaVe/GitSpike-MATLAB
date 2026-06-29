%% Script to generate a dataset in labeled line population 
% Date:June 2026
% Author : Laure WOLFF

function CellMatrix = generate_and_plot_raster_ll(num_stimuli, num_repetitions, ...
    num_indi, num_neurons, t1, t2, base_rate, refrac, showing, plotting, other_figs)
    
    % Input parameter consistency check
    if num_indi > num_neurons
        error('Error: The number of the individuel coding neuron cannot exceed the total num_neurons!');
    end
    
    num_trials = num_stimuli * num_repetitions;
    CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);

    %% 1. Generation of INDIVIDUAL coding neurons
    % Define a preference matrix to map individual neuron activations to stimuli
    % Rows = Neurons, Columns = Stimuli. 
    % A value of 1 means the neuron responds to that specific stimulus.
    if num_indi == 4 && num_stimuli == 4 
        pref_matrix = [
            1, 1, 0, 0;  % Neuron 1: Sensitive to Feature A1 (e.g., Color White) -> fires for S1, S2
            0, 0, 1, 1;  % Neuron 2: Sensitive to Feature A2 (e.g., Color Red)   -> fires for S3, S4
            1, 0, 1, 0;  % Neuron 3: Sensitive to Feature B1 (e.g., Car)         -> fires for S1, S3
            0, 1, 0, 1   % Neuron 4: Sensitive to Feature B2 (e.g., Ship)        -> fires for S2, S4
        ];
    else
        % Fallback for safety if you change sizes later, but randi is not ideal for the benchmark
        pref_matrix = randi([0 1], num_indi, num_stimuli);
    end

    if showing == true
        fprintf('\n--- Labeled Line Preference Matrix (Neurons x Stimuli) ---\n');
        disp(pref_matrix);
        fprintf('----------------------------------------------------------\n\n');
    end
    for c_idx = 1:num_indi
        for st = 1:num_stimuli
            
            % If neuron 'c_idx' prefers stimulus 'st', increase its firing rate (Rate Coding)
            % Otherwise, it behaves like a non-coding background noise neuron (base_rate)
            if pref_matrix(mod(c_idx-1, size(pref_matrix,1))+1, st) == 1
                local_rate = base_rate * 3.0; % High activity (selective coding)
            else
                local_rate = base_rate * 1.0; % Baseline noise level
            end
            
            for rp = 1:num_repetitions
                approx_spikes = round((t2 - t1) * local_rate * 3) + 10;
                uniform_samples = rand(1, approx_spikes);
                intervals = refrac - log(1 - uniform_samples) / local_rate;
                spikes = cumsum(intervals);
                CellMatrix{c_idx, st, rp} = spikes(spikes >= t1 & spikes <= t2);
            end
        end
    end
    
    %% 3. Generation of NON-CODING neurons (Channels: num_indi + 1 to num_neurons)
    for st = 1:num_stimuli
        for rp = 1:num_repetitions
            for c_idx = num_indi + 1:num_neurons
                approx_spikes = round((t2 - t1) * base_rate * 1.0) + 10;
                uniform_samples = rand(1, approx_spikes);
                intervals = refrac - log(1 - uniform_samples) / base_rate;
                spikes_noise = cumsum(intervals);
                CellMatrix{c_idx, st, rp} = spikes_noise(spikes_noise >= t1 & spikes_noise <= t2);
            end
        end
    end
    
    %% =========================================================================
    %% PLOTTING SECTION (Layout matching the 2018 Paper Style)
    %% =========================================================================
    if plotting == true
        if other_figs == true
            figure('Name', 'Individual Neuronal Raster Plots (LL Mode)', 'Color', 'w', 'Position', [100, 100, 1400, 900]);
    
    % Automatically calculate the grid size (e.g., 4 rows x 5 columns for 20 neurons)
    cols = ceil(sqrt(num_neurons * 1.25)); 
    rows = ceil(num_neurons / cols);
    
    % Total number of trials (lines on the Y-axis for each subplot)
    num_trials = num_stimuli * num_repetitions;
    
    %% Loop over each individual neuron
    for n = 1:num_neurons
        subplot(rows, cols, n);
        hold on;
        
        % Background colors to visually separate stimuli blocks on the raster
        for st = 1:num_stimuli
            % Define Y-range for the current stimulus block
            y_start = (st - 1) * num_repetitions + 0.5;
            y_end = st * num_repetitions + 0.5;
            
            % Alternating background colors (light gray and white) for contrast
            if mod(st, 2) == 1
                fill([t1 t2 t2 t1], [y_start y_start y_end y_end], [0.96 0.96 0.96], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
            end
        end
        
        % Plot the spikes for this specific neuron
        t_counter = 1;
        for st = 1:num_stimuli
            for rp = 1:num_repetitions
                spikes = CellMatrix{n, st, rp};
                
                if ~isempty(spikes)
                    % Draw a small vertical line for each spike
                    for sp = 1:length(spikes)
                        line([spikes(sp), spikes(sp)], [t_counter - 0.4, t_counter + 0.4], ...
                             'Color', 'k', 'LineWidth', 1.0);
                    end
                end
                t_counter = t_counter + 1;
            end
        end
        
        % Draw horizontal lines to cleanly separate the stimuli blocks
        for st_sep = 1:(num_stimuli-1)
            sep_line = st_sep * num_repetitions + 0.5;
            line([t1, t2], [sep_line, sep_line], 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0, 'LineStyle', '--');
        end
        
        %% Layout and Axis formatting for each subplot
        box on; 
        xlim([t1, t2]); 
        ylim([0.5, num_trials + 0.5]);
        set(gca, 'YDir', 'reverse'); % Traditional raster orientation (Trial 1 at the top)
        
        % Tick labels configuration
        if num_stimuli <= 5
            % Place a tick in the middle of each stimulus block
            y_ticks = (0.5 : num_repetitions : num_trials) + (num_repetitions / 2);
            y_labels = cell(1, num_stimuli);
            for s = 1:num_stimuli
                y_labels{s} = sprintf('S%d', s);
            end
            set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels, 'FontSize', 7);
        else
            set(gca, 'YTick', [1, num_trials], 'YTickLabel', {'1', num2str(num_trials)}, 'FontSize', 7);
        end
        
        % Title customized by neuron type (Red for Coding, Blue for Noise)
        if n <= num_indi
            title(sprintf('Neuron %d (Coding)', n), 'Color', [0.85 0 0], 'FontSize', 9, 'FontWeight', 'bold');
        else
            title(sprintf('Neuron %d (Noise)', n), 'Color', [0 0.45 0.74], 'FontSize', 9, 'FontWeight', 'normal');
        end
        
        % Add axis labels only to the outer subplots to avoid text clutter
        if mod(n-1, cols) == 0
            ylabel('Stimuli / Trials', 'FontSize', 7);
        end
        if n > (num_neurons - cols)
            xlabel('Time', 'FontSize', 7);
        end
    end
    
    hold off;
        end

        %% Multi-Subplot layout per Trial (Pooled lines for C, NC and All at the bottom)
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
                    if c_idx <=  num_indi
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
                    y_pos = (num_neurons - c_idx + 1) + 4; % Inversion pour avoir le neurone 1 en haut
                    
                    if c_idx <=  num_indi
                        current_color = [1 0 0]; % Red
                        len = length(spikes);
                        if len > 0
                            all_indi_spikes(ptr_i : ptr_i + len - 1) = spikes;
                            ptr_i = ptr_i + len;
                        end
                    else
                        current_color = [0.0000 0.4470 0.7410]; % Blue
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
                
                % Data cleanup and chronological sorting
                all_indi_spikes  = unique(sort(all_indi_spikes(1:ptr_i-1)));
                all_noise_spikes = unique(sort(all_noise_spikes(1:ptr_nc-1)));
                all_total_spikes = unique(sort([all_indi_spikes, all_noise_spikes]));
                
                % Plotting the integrated summation tracks at the bottom
                for sp = 1:length(all_total_spikes)
                    line([all_total_spikes(sp), all_total_spikes(sp)], [1 - 0.4, 1 + 0.4], 'Color', [0 0 0], 'LineWidth', 1.3);
                end
                for sp = 1:length(all_noise_spikes)
                    line([all_noise_spikes(sp), all_noise_spikes(sp)], [2 - 0.3, 2 + 0.3], 'Color', [0.0000 0.4470 0.7410], 'LineWidth', 1.0);
                end
                for sp = 1:length(all_indi_spikes)
                    line([all_indi_spikes(sp), all_indi_spikes(sp)], [3 - 0.3, 3 + 0.3], 'Color', [1 0 0], 'LineWidth', 1.5);
                end
               
                
                % Visual dashed separators for the subpopulation blocks
                line([t1, t2], [4.5, 4.5], 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
                line([t1, t2], [num_neurons + 4.5 -  num_indi, num_neurons + 4.5 - num_indi], 'Color', [0.7 0.7 0.7], 'LineStyle', ':');
                
                box on; xlim([t1, t2]); ylim([0.5, num_neurons + 5.5]);
                
                % Define ticks and labels dynamically (Maintained robustness for num_coll=0 or num_indi=0)
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