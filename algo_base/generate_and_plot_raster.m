%% Script to generate a dataset in summed populatio like tell in the 2018 paper
% Date: May-June 2026
% Author : Laure WOLFF

function CellMatrix = generate_and_plot_raster(num_stimuli, ...
    num_repetitions, num_coding_neurons, num_neurons, t1, t2, base_rate, refrac, plotting,other_figs)
    num_trials = num_stimuli * num_repetitions;
    CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);
        
    %% 1. Generation of coding neurons (1 to c)
    for st = 1:num_stimuli
        % The global pooled spike train for stimulus 'st' must have a total rate 
        % equal to (c * base_rate) so that after distribution, each neuron has a rate of M.
        pooled_rate = num_coding_neurons * base_rate;
        
        % Generate a homogeneous Poisson process with an absolute refractory period
        approx_spikes = round((t2 - t1) * pooled_rate * 3) + 10;
        uniform_samples = rand(1, approx_spikes);
        intervals = refrac - log(1 - uniform_samples) / pooled_rate;
        spikes_pooled = cumsum(intervals);
        spikes_pooled = spikes_pooled(spikes_pooled >= t1 & spikes_pooled <= t2);
        
        num_spikes = length(spikes_pooled);
        
        % For each repetion, distribution of the spikes between each
        % coding-neurons trains
        for rp = 1:num_repetitions
            if num_spikes > 0
                % Random but fair distribution of spike indices.
                % (Each neuron receives an exclusive subset of spikes for each trial)
                shuffled_indices = randperm(num_spikes);
                
                for nc = 1:num_coding_neurons
                    % Select the spikes assigned to neuron 'nc' for this trial 'rp'
                    idx_assigned = shuffled_indices(nc:num_coding_neurons:end);
                    CellMatrix{nc, st, rp} = sort(spikes_pooled(idx_assigned));
                end
            else
                for nc = 1:num_coding_neurons
                    CellMatrix{nc, st, rp} = [];
                end
            end
        end
    end
    
    %% Generation of non coding-neurons (c+1 to N)
    % For these, each trial of every stimulus is an independent Poisson process with rate M.
    for st = 1:num_stimuli
        for rp = 1:num_repetitions
            for nc = (num_coding_neurons + 1):num_neurons
                approx_spikes = round((t2 - t1) * base_rate * 3) + 10;
                uniform_samples = rand(1, approx_spikes);
                intervals = refrac - log(1 - uniform_samples) / base_rate;
                
                spikes_noise = cumsum(intervals);
                CellMatrix{nc, st, rp} = spikes_noise(spikes_noise >= t1 & spikes_noise <= t2);
            end
        end
    end
    
    %% Plotting 
    if plotting == true
        if other_figs == true
            color_coding  = [0.85 0.33 0.1];  % Red/Orange for coding neurons (C)
            color_noise   = [0.0 0.45 0.74];  % Blue for non-coding neurons (NC)

            plot_types = {'CODING', 'NON-CODING', 'FULL'};
            titles = { ...
                sprintf('1. Coding Subpopulation Only (Neurons 1 to %d)', num_coding_neurons), ...
                sprintf('2. Non-Coding Subpopulation Only (Neurons %d to %d)', num_coding_neurons + 1, num_neurons), ...
                sprintf('3. Full Population (All %d Neurons)', num_neurons) ...
            };

            % Generate consistent labels for the trials
            trial_labels = cell(1, num_trials);
            counter = 1;
            for st = 1:num_stimuli
                for rp = 1:num_repetitions
                    trial_labels{counter} = sprintf('S%d-R%d', st, rp);
                    counter = counter + 1;
                end
            end

            % Plotting
            figure('Name', 'Subpopulations spike train raster plots', 'Position', [200, 50, 950, 950]);

            for f = 1:3
                subplot(3, 1, f);
                hold on;

                for t_idx = 1:num_trials
                    st = floor((t_idx-1)/num_repetitions) + 1;
                    rp = mod((t_idx-1), num_repetitions) + 1;

                    if strcmp(plot_types{f}, 'CODING')
                        neurons_to_plot = 1:num_coding_neurons;
                    elseif strcmp(plot_types{f}, 'NON-CODING')
                        neurons_to_plot = (num_coding_neurons + 1):num_neurons;
                    else
                        neurons_to_plot = 1:num_neurons;
                    end

                    for nc = neurons_to_plot
                        spikes = CellMatrix{nc, st, rp};
                        if ~isempty(spikes)
                            if nc <= num_coding_neurons
                                current_color = color_coding; line_width = 1.5;
                            else
                                current_color = color_noise; line_width = 1.0;
                            end
                            for sp = 1:length(spikes)
                                % Standard trial line plotting
                                line([spikes(sp), spikes(sp)], [t_idx - 0.35, t_idx + 0.35], ...
                                     'Color', current_color, 'LineWidth', line_width);
                            end
                        end
                    end
                    line([t1, t2], [t_idx, t_idx], 'Color', [0.94 0.94 0.94], 'LineWidth', 0.5);
                end

                for st_sep = 1:(num_stimuli-1)
                    sep_line = st_sep * num_repetitions + 0.5;
                    line([t1, t2], [sep_line, sep_line], 'Color', [0.1 0.1 0.1], 'LineWidth', 1.3, 'LineStyle', '--');
                end

                box on; grid on;
                set(gca, 'XGrid', 'on', 'YGrid', 'off', 'YDir', 'reverse'); % Invert the Y-axis to have S1-R1 at the top
                xlim([t1, t2]); ylim([0.5, num_trials + 0.5]);

                set(gca, 'YTick', 1:num_trials, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none', 'FontSize', 8);
                ylabel('Trials (Stimuli & Reps)', 'FontSize', 10, 'FontWeight', 'bold');
                title(titles{f}, 'FontSize', 11, 'FontWeight', 'bold');

                if f == 3
                    xlabel('Time (au)', 'FontSize', 11, 'FontWeight', 'bold');
                end
                hold off;
            end

            % trial_labels = cell(1, num_trials);
            % counter = 1;
            % for st = 1:num_stimuli
            %     for rp = 1:num_repetitions
            %         trial_labels{counter} = sprintf('S%d-R%d', st, rp);
            %         counter = counter + 1;
            %     end
            % end
            % 
            % % Plotting the global raster plot
            % figure('Name', 'Global dataset raster plot', 'Position', [100, 100, 950, 650]);
            % hold on;
            % 
            % color_noise  = [0.0 0.45 0.74];  % Blue for the non-coding neurons (NC)
            % color_coding = [0.85 0.33 0.1];  % Red for the coding neurons (C)
            % 
            % for t_idx = 1:num_trials
            %     st = floor((t_idx-1)/num_repetitions) + 1;
            %     rp = mod((t_idx-1), num_repetitions) + 1;
            % 
            %     for nc = 1:num_neurons
            %         spikes = CellMatrix{nc, st, rp};
            %         if ~isempty(spikes)
            %             if nc <= num_coding_neurons
            %                 current_color = color_coding; line_width = 1.5;
            %             else
            %                 current_color = color_noise; line_width = 1.0;
            %             end
            %             for sp = 1:length(spikes)
            %                 % Plots each trial row
            %                 line([spikes(sp), spikes(sp)], [t_idx - 0.4, t_idx + 0.4], ...
            %                      'Color', current_color, 'LineWidth', line_width);
            %             end
            %         end
            %     end
            %     line([t1, t2], [t_idx, t_idx], 'Color', [0.94 0.94 0.94], 'LineWidth', 0.5, 'LineStyle', '-');
            % end
            % 
            % set(gca, 'YDir', 'reverse'); % Invert the Y-axis to have S1-R1 at the top
            % set(gca, 'XGrid', 'on', 'YGrid', 'off');
            % xlim([t1, t2]); 
            % ylim([0.5, num_trials + 0.5]);
            % 
            % set(gca, 'YTick', 1:num_trials, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none', 'FontSize', 9);
            % xlabel('Time', 'FontSize', 11, 'FontWeight', 'bold');
            % ylabel('Trials (Stimuli & Repetitions)', 'FontSize', 11, 'FontWeight', 'bold');
            % title('Global spike train raster plot (Red: Coding | Blue: Noise)', 'FontSize', 13, 'FontWeight', 'bold');
            % box on; 
            % grid on;
            % 
            % for st = 1:(num_stimuli-1)
            %     sep_line = st * num_repetitions + 0.5;
            %     line([t1, t2], [sep_line, sep_line], 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2, 'LineStyle', '--');
            % end
            shg;
        end

        % Creation of the Figure 1 of the 2018's paper 
        figure('Name', 'Subpopulations pooled spike train raster plots');
        idx_subplot = 0; 
        for st_select = 1:num_stimuli
            for rp_select = 1:num_repetitions
                idx_subplot = idx_subplot + 1;
                subplot(num_trials, 1, idx_subplot);
                hold on;
                max_coding_size = 0;
                max_noise_size = 0;
                for nc = 1:num_neurons
                    sz = length(CellMatrix{nc, st_select, rp_select});
                    if nc <= num_coding_neurons
                        max_coding_size = max_coding_size + sz;
                    else
                        max_noise_size = max_noise_size + sz;
                    end
                end

                all_coding_spikes = zeros(1, max_coding_size);
                all_noise_spikes = zeros(1, max_noise_size);
                ptr_c = 1;
                ptr_nc = 1;

                % plotting the non-coding and coding neurons
                for nc = 1:num_neurons
                    spikes = CellMatrix{nc, st_select, rp_select};
                    y_pos = num_neurons - nc + 4; 

                    if nc <= num_coding_neurons
                        current_color = [0.85 0.33 0.1]; % Red for the coding neurons
                        len = length(spikes);
                        if len > 0
                            all_coding_spikes(ptr_c : ptr_c + len - 1) = spikes;
                            ptr_c = ptr_c + len;
                        end
                    else
                        current_color = [0.0 0.45 0.74]; % Blue for non-coding neurons 
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

                % Several lines to clearify the plot
                line([t1, t2], [num_neurons+3.5-num_coding_neurons, num_neurons+3.5-num_coding_neurons], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 0.8);
                line([t1, t2], [3.5, 3.5], 'Color', [0 0 0], 'LineWidth', 1.2);

                % adding all spike train of each neuron to plot the line c, cn  and all
                all_coding_spikes = unique(sort(all_coding_spikes(1:ptr_c-1)));
                all_noise_spikes = unique(sort(all_noise_spikes(1:ptr_nc-1)));
                all_total_spikes = unique(sort([all_coding_spikes, all_noise_spikes]));

                % Line C 
                for sp = 1:length(all_coding_spikes)
                    line([all_coding_spikes(sp), all_coding_spikes(sp)], [3 - 0.3, 3 + 0.3], ...
                         'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
                end

                % Line NC
                for sp = 1:length(all_noise_spikes)
                    line([all_noise_spikes(sp), all_noise_spikes(sp)], [2 - 0.3, 2 + 0.3], ...
                         'Color', [0.0 0.45 0.74], 'LineWidth', 1.0);
                end

                % Line All
                for sp = 1:length(all_total_spikes)
                    line([all_total_spikes(sp), all_total_spikes(sp)], [1 - 0.4, 1 + 0.4], ...
                         'Color', [0 0 0], 'LineWidth', 1.2);
                end

                box on;
                xlim([t1, t2]); 
                ylim([0.5, num_neurons + 4.5]);
                y_ticks = [1, 2, 3, 4, num_neurons+4-num_coding_neurons, num_neurons+3];
                y_labels = {'All', 'NC', 'C', num2str(num_neurons), num2str(num_coding_neurons), '1'};
                set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels, 'FontSize', 7);
                title(sprintf('trial : S%d-R%d', st_select, rp_select), ...
                    'FontSize', 8, 'FontWeight', 'bold');
                set(gca, 'XTickLabel', []);
                ylabel('Spike trains', 'FontSize', 8);

            end
        end
        xlabel('Time (au)', 'FontSize', 8);
    end

        
        
        