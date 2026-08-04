%% Wrapper MATLAB pour appeler le moteur C-MEX et afficher les Rasters
% Date: Mai-Juin 2026
% Author : Laure WOLFF

function CellMatrix = generate_and_plot_raster_mex_m(num_stimuli, num_repetitions, num_coding_neurons, num_neurons, t1, t2, base_rate, refrac, plotting, other_figs)
    
    %% Call of the MEX file
    CellMatrix = generate_and_plot_raster_mex(num_stimuli, num_repetitions, num_coding_neurons, num_neurons, t1, t2, base_rate, refrac);
    
    %% Plotting
    if plotting == true
        num_trials = num_stimuli * num_repetitions;
        if other_figs == true
            color_coding  = [0.85 0.33 0.1];  
            color_noise   = [0.0 0.45 0.74];  

            plot_types = {'CODING', 'NON-CODING', 'FULL'};
            titles = { ...
                sprintf('1. Coding Subpopulation Only (Neurons 1 to %d)', num_coding_neurons), ...
                sprintf('2. Non-Coding Subpopulation Only (Neurons %d to %d)', num_coding_neurons + 1, num_neurons), ...
                sprintf('3. Full Population (All %d Neurons)', num_neurons) ...
            };

            trial_labels = cell(1, num_trials);
            counter = 1;
            for st = 1:num_stimuli
                for rp = 1:num_repetitions
                    trial_labels{counter} = sprintf('S%d-R%d', st, rp);
                    counter = counter + 1;
                end
            end

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
                set(gca, 'XGrid', 'on', 'YGrid', 'off', 'YDir', 'reverse'); 
                xlim([t1, t2]); ylim([0.5, num_trials + 0.5]);
                set(gca, 'YTick', 1:num_trials, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none', 'FontSize', 8);
                ylabel('Trials (Stimuli & Reps)', 'FontSize', 10, 'FontWeight', 'bold');
                title(titles{f}, 'FontSize', 11, 'FontWeight', 'bold');
                if f == 3, xlabel('Time (au)', 'FontSize', 11, 'FontWeight', 'bold'); end
                hold off;
            end
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
                max_coding_size = 0; max_noise_size = 0;
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
                ptr_c = 1; ptr_nc = 1;

                for nc = 1:num_neurons
                    spikes = CellMatrix{nc, st_select, rp_select};
                    y_pos = num_neurons - nc + 4; 

                    if nc <= num_coding_neurons
                        current_color = [0.85 0.33 0.1]; 
                        len = length(spikes);
                        if len > 0
                            all_coding_spikes(ptr_c : ptr_c + len - 1) = spikes;
                            ptr_c = ptr_c + len;
                        end
                    else
                        current_color = [0.0 0.45 0.74]; 
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

                line([t1, t2], [num_neurons+3.5-num_coding_neurons, num_neurons+3.5-num_coding_neurons], 'Color', [0.6 0.6 0.6], 'LineStyle', '--', 'LineWidth', 0.8);
                line([t1, t2], [3.5, 3.5], 'Color', [0 0 0], 'LineWidth', 1.2);

                all_coding_spikes = unique(sort(all_coding_spikes(1:ptr_c-1)));
                all_noise_spikes = unique(sort(all_noise_spikes(1:ptr_nc-1)));
                all_total_spikes = unique(sort([all_coding_spikes, all_noise_spikes]));

                for sp = 1:length(all_coding_spikes)
                    line([all_coding_spikes(sp), all_coding_spikes(sp)], [3 - 0.3, 3 + 0.3], 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
                end
                for sp = 1:length(all_noise_spikes)
                    line([all_noise_spikes(sp), all_noise_spikes(sp)], [2 - 0.3, 2 + 0.3], 'Color', [0.0 0.45 0.74], 'LineWidth', 1.0);
                end
                for sp = 1:length(all_total_spikes)
                    line([all_total_spikes(sp), all_total_spikes(sp)], [1 - 0.4, 1 + 0.4], 'Color', [0 0 0], 'LineWidth', 1.2);
                end

                box on; xlim([t1, t2]); ylim([0.5, num_neurons + 4.5]);
                y_ticks = [1, 2, 3, 4, num_neurons+4-num_coding_neurons, num_neurons+3];
                y_labels = {'All', 'NC', 'C', num2str(num_neurons), num2str(num_coding_neurons), '1'};
                set(gca, 'YTick', y_ticks, 'YTickLabel', y_labels, 'FontSize', 7);
                title(sprintf('trial : S%d-R%d', st_select, rp_select), 'FontSize', 8, 'FontWeight', 'bold');
                set(gca, 'XTickLabel', []); ylabel('Spike trains', 'FontSize', 8);
            end
        end
        xlabel('Time (au)', 'FontSize', 8);
    end
end