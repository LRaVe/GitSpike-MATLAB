%% Script to generate a dataset in summed populatio like tell in the 2018 paper
% Date: May-June 2026
% Author : Laure WOLFF

function CellMatrix = generate_and_plot_raster(num_stimuli, num_repetitions, num_coding_neurons, num_neurons, t1, t2, base_rate, refrac, plotting)
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
        color_coding  = [0.85 0.33 0.1];  % Orange/Rouge pour les codants
        color_noise   = [0.0 0.45 0.74];  % Bleu pour le bruit
    
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
    
        figure('Name', 'Subpopulations Spike Train Raster Plots (Satuvuori SP Method)', 'Position', [200, 50, 950, 950]);
    
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
            set(gca, 'XGrid', 'on', 'YGrid', 'off');
            xlim([t1, t2]); ylim([0.5, num_trials + 0.5]);
            
            set(gca, 'YTick', 1:num_trials, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none', 'FontSize', 8);
            ylabel('Trials (Stimuli & Reps)', 'FontSize', 10, 'FontWeight', 'bold');
            title(titles{f}, 'FontSize', 11, 'FontWeight', 'bold');
            
            if f == 3
                xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
            end
            hold off;
        end
        shg;
    end
end