%% Script to generate a dataset in summed population like told in the 2018 paper
% Special Case: Individual vs. Collective (Coll) vs. Non-Coding (NC)
% Date: May-June 2026
% Author : Laure WOLFF

function CellMatrix = generate_and_plot_raster(num_stimuli, num_repetitions, ...
    num_indi, num_coll, num_neurons, t1, t2, base_rate, refrac, plotting, other_figs)
    
    % Input parameter consistency check
    if (num_indi + num_coll) > num_neurons
        error('Error: The sum of num_indi and num_coll cannot exceed the total num_neurons!');
    end
    
    num_trials = num_stimuli * num_repetitions;
    CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);
    
    %% 1. Generation of INDIVIDUAL coding neurons (Channels: 1 to num_indi)
    % Each individual neuron carries its own exclusive discriminant information
    for nc = 1:num_indi
        for st = 1:num_stimuli
            % Shift the rate slightly depending on the stimulus to create distinct coding properties
            local_rate = base_rate * (1 + 0.25 * (st - 1)); 
            
            for rp = 1:num_repetitions
                approx_spikes = round((t2 - t1) * local_rate * 3) + 10;
                uniform_samples = rand(1, approx_spikes);
                intervals = refrac - log(1 - uniform_samples) / local_rate;
                spikes = cumsum(intervals);
                CellMatrix{nc, st, rp} = spikes(spikes >= t1 & spikes <= t2);
            end
        end
    end
    
    %% 2. Generation of COLLECTIVE coding neurons (Channels: num_indi + 1 to num_indi + num_coll)
    % Core concept of the 2018 paper: Taken individually, each neuron looks like noise.
    % However, the SUM (union) of their spikes reconstructs a unique pattern per stimulus.
    for st = 1:num_stimuli
        % The collective pool holds the structural pattern of the stimulus
        pooled_rate = num_coll * base_rate * 1.1; 
        
        approx_spikes = round((t2 - t1) * pooled_rate * 3) + 10;
        uniform_samples = rand(1, approx_spikes);
        intervals = refrac - log(1 - uniform_samples) / pooled_rate;
        spikes_pooled = cumsum(intervals);
        spikes_pooled = spikes_pooled(spikes_pooled >= t1 & spikes_pooled <= t2);
        num_spikes = length(spikes_pooled);
        
        for rp = 1:num_repetitions
            if num_spikes > 0
                shuffled_indices = randperm(num_spikes);
                for c_idx = 1:num_coll
                    nc = num_indi + c_idx; % Channels mapped right after individual neurons
                    idx_assigned = shuffled_indices(c_idx:num_coll:end);
                    CellMatrix{nc, st, rp} = sort(spikes_pooled(idx_assigned));
                end
            else
                for c_idx = 1:num_coll
                    CellMatrix{num_indi + c_idx, st, rp} = [];
                end
            end
        end
    end
    
    %% 3. Generation of NON-CODING neurons (NC Channels: from num_indi + num_coll + 1 to num_neurons)
    % Pure independent noise background, without any stimulus-driven structure.
    for st = 1:num_stimuli
        for rp = 1:num_repetitions
            for nc = (num_indi + num_coll + 1):num_neurons
                approx_spikes = round((t2 - t1) * base_rate * 1.0) + 10;
                uniform_samples = rand(1, approx_spikes);
                intervals = refrac - log(1 - uniform_samples) / base_rate;
                spikes_noise = cumsum(intervals);
                CellMatrix{nc, st, rp} = spikes_noise(spikes_noise >= t1 & spikes_noise <= t2);
            end
        end
    end
    
    %% =========================================================================
    %% PLOTTING SECTION (Layout matching the 2018 Paper Style)
    %% =========================================================================
    if plotting == true
        if other_figs == true
            % Distinct color scheme to visually identify the 3 subgroups from the paper
            color_indi = [1 0 0];  % Red for Individual (Indi)
            color_coll = [0.8500 0.3250 0.0980];  % Orange/Red for Collective (Coll)
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

                for nc = 1:num_neurons
                    spikes = CellMatrix{nc, st, rp};
                    if ~isempty(spikes)
                        % Assign color boundaries based on custom input dimensions
                        if nc <= num_indi
                            current_color = color_indi; line_width = 1.4;
                        elseif nc <= (num_indi + num_coll)
                            current_color = color_coll; line_width = 1.4;
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

            % Horizontal boundary lines separating different stimuli blocks
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
            title('Artificial Dataset (Red: Indi | Orange: Coll | Blue: NC Noise)', 'FontSize', 12, 'FontWeight', 'bold');
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
                
                max_indi_size = 0; max_coll_size = 0; max_noise_size = 0;
                for nc = 1:num_neurons
                    sz = length(CellMatrix{nc, st_select, rp_select});
                    if nc <= num_indi, max_indi_size = max_indi_size + sz;
                        elseif nc <= (num_indi + num_coll), max_coll_size = max_coll_size + sz;
                        else, max_noise_size = max_noise_size + sz; 
                    end
                end

                all_indi_spikes = zeros(1, max_indi_size);
                all_coll_spikes = zeros(1, max_coll_size);
                all_noise_spikes = zeros(1, max_noise_size);
                
                ptr_i = 1; ptr_co = 1; ptr_nc = 1;

                for nc = 1:num_neurons
                    spikes = CellMatrix{nc, st_select, rp_select};
                    y_pos = num_neurons - nc + 5; 

                    if nc <= num_indi
                        current_color = [1 0 0];
                        len = length(spikes);
                        if len > 0
                            all_indi_spikes(ptr_i : ptr_i + len - 1) = spikes;
                            ptr_i = ptr_i + len;
                        end
                    elseif nc <= (num_indi + num_coll)
                        current_color = [0.8500 0.3250 0.0980];
                        len = length(spikes);
                        if len > 0
                            all_coll_spikes(ptr_co : ptr_co + len - 1) = spikes;
                            ptr_co = ptr_co + len;
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

                % Data cleanup and chronological sorting of unique accumulated vector points
                all_indi_spikes  = unique(sort(all_indi_spikes(1:ptr_i-1)));
                all_coll_spikes  = unique(sort(all_coll_spikes(1:ptr_co-1)));
                all_noise_spikes = unique(sort(all_noise_spikes(1:ptr_nc-1)));
                all_total_spikes = unique(sort([all_indi_spikes, all_coll_spikes, all_noise_spikes]));

                % Plotting the integrated summation tracks at the bottom of the axis
                for sp = 1:length(all_indi_spikes)
                    line([all_indi_spikes(sp), all_indi_spikes(sp)], [4 - 0.3, 4 + 0.3], 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.2);
                end
                for sp = 1:length(all_coll_spikes)
                    line([all_coll_spikes(sp), all_coll_spikes(sp)], [3 - 0.3, 3 + 0.3], 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 1.2);
                end
                for sp = 1:length(all_noise_spikes)
                    line([all_noise_spikes(sp), all_noise_spikes(sp)], [2 - 0.3, 2 + 0.3], 'Color', [0.0000 0.4470 0.7410], 'LineWidth', 1.0);
                end
                for sp = 1:length(all_total_spikes)
                    line([all_total_spikes(sp), all_total_spikes(sp)], [1 - 0.4, 1 + 0.4], 'Color', [0 0 0], 'LineWidth', 1.3);
                end

                % Visual dashed separators for the subpopulation blocks
                line([t1, t2], [4.5, 4.5], 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2);
                line([t1, t2], [num_neurons+4.5-(num_indi+num_coll), num_neurons+4.5-(num_indi+num_coll)], 'Color', [0.7 0.7 0.7], 'LineStyle', ':');
                line([t1, t2], [num_neurons+4.5-num_indi, num_neurons+4.5-num_indi], 'Color', [0.7 0.7 0.7], 'LineStyle', ':');

                box on; xlim([t1, t2]); ylim([0.5, num_neurons + 5.5]);
                
                % Define ticks and labels dynamically
                y_ticks = [1, 2, 3, 4, 5, num_neurons+4-num_indi, num_neurons+4];
                y_labels = {'Total', '\Sigma NC', '\Sigma Coll', '\Sigma Indi', num2str(num_neurons), num2str(num_indi+1), '1'};
                
                % Fix for num_indi = 0 or specific cases creating overlapping ticks
                [y_ticks, unique_idx] = unique(y_ticks);
                y_labels = y_labels(unique_idx);
                
                % Apply the filtered and strictly increasing ticks to the axes
                set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels, 'FontSize', 7);
                
                title(sprintf('Trial : S%d-R%d', st_select, rp_select), 'FontSize', 8, 'FontWeight', 'bold');
                set(gca, 'XTickLabel', []);
            end
        end
        xlabel('Time (au)', 'FontSize', 9, 'FontWeight', 'bold');
    end
end