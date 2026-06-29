%% Script to generate a dataset in summed population like told in the 2018 paper
% Special Case: Collective (Coll) vs. Individual (Indi) vs. Non-Coding (NC)
% Date:June 2026
% Author : Laure WOLFF

function CellMatrix = generate_and_plot_raster_ll(num_stimuli, num_repetitions, ...
    num_indi, num_neurons, t1, t2, base_rate, refrac, plotting, other_figs)
    
    % Input parameter consistency check
    if num_indi > num_neurons
        error('Error: The number of the individuel coding neuron cannot exceed the total num_neurons!');
    end
    
    num_trials = num_stimuli * num_repetitions;
    CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);

    %% 1. Generation of INDIVIDUAL coding neuronz
    for c_idx = 1:num_indi

        for st = 1:num_stimuli
            % On crée un profil de fréquence contrasté par stimulus
            % Stimulus 1 -> base_rate * 0.75 (Fréquence basse)
            % Stimulus 2 -> base_rate * 1.75 (Fréquence haute)
            local_rate = base_rate * (0.75 + 1.0 * (st - 1)); 
            
            for rp = 1:num_repetitions
                approx_spikes = round((t2 - t1) * local_rate * 3) + 10;
                uniform_samples = rand(1, approx_spikes);
                intervals = refrac - log(1 - uniform_samples) / local_rate;
                spikes = cumsum(intervals);
                CellMatrix{c_idx, st, rp} = spikes(spikes >= t1 & spikes <= t2);
            end
        end
    end
    
    %% 3. Generation of NON-CODING neurons (Channels: num_coll + num_indi + 1 to num_neurons)
    for st = 1:num_stimuli
        for rp = 1:num_repetitions
            for c_idx = (num_coll + num_indi + 1):num_neurons
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
            % Distinct color scheme to visually identify the 2 subgroups
            color_indi = [1 0 0];                 % Red for Individual (Indi)
            color_noise = [0.0000 0.4470 0.7410]; % Blue for Non-Coding (NC)
            trial_labels = cell(1, num_trials);
            counter = 1;
            for st = 1:num_stimuli
                for rp = 1:num_repetitions
                    trial_labels{counter} = sprintf('S%d-R%d', st, rp);
                    counter = counter + 1;
                end
            end
            figure('Name', ' Global Raster Plot', 'Color', 'w', 'Position', [100, 100, 950, 700]);
            hold on;
            
            for t_idx = 1:num_trials
                st = floor((t_idx-1)/num_repetitions) + 1;
                rp = mod((t_idx-1), num_repetitions) + 1;
                for c_idx = 1:num_neurons
                    spikes = CellMatrix{c_idx, st, rp};
                    if ~isempty(spikes)
                        if c_idx <= num_indi
                            current_color = color_indi; line_width = 1.5;
                        else
                            current_color = color_noise; line_width = 1.0;
                        end
                        
                        for sp = 1:length(spikes)
                            line([spikes(sp), spikes(sp)], [t_idx - 0.35, t_idx + 0.35], ...
                                 'Color', current_color, 'LineWidth', line_width);
                        end
                    end
                end
                line([t1, t2], [t_idx, t_idx], 'Color', [0.95 0.95 0.95], 'LineWidth', 0.5);
            end
            for st_sep = 1:(num_stimuli-1)
                sep_line = st_sep * num_repetitions + 0.5;
                line([t1, t2], [sep_line, sep_line], 'Color', [0.2 0.2 0.2], 'LineWidth', 1.5, 'LineStyle', '--');
            end
            box on; grid on;
            set(gca, 'XGrid', 'on', 'YGrid', 'off', 'YDir', 'reverse');
            xlim([t1, t2]); ylim([0.5, num_trials + 0.5]);
            set(gca, 'YTick', 1:num_trials, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none', 'FontSize', 9);
            
            xlabel('Time (au)', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('Trials (Stimuli / Repetitions)', 'FontSize', 11, 'FontWeight', 'bold');
            title('Artificial Dataset ( Red: Indi | Blue: NC Noise)', 'FontSize', 12, 'FontWeight', 'bold');
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