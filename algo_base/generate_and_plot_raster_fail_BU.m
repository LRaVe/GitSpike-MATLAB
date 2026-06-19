%% Script to generate a dataset in summed population like told in the 2018 paper
% Special Case: Collective (Coll) vs. Individual (Indi) vs. Non-Coding (NC)
% Date: May-June 2026
% Author : Laure WOLFF
function CellMatrix = generate_and_plot_raster_fail_BU(num_stimuli, num_repetitions, ...
    num_indi, num_coll, num_neurons, t1, t2, base_rate, refrac, plotting, other_figs)
    
    % Input parameter consistency check
    if (num_indi + num_coll) > num_neurons
        error('Error: The sum of num_indi and num_coll cannot exceed the total num_neurons!');
    end
    
    num_trials = num_stimuli * num_repetitions;
    CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);
    
    %% 1. Generation of COLLECTIVE coding neurons (Channels: 1 to num_coll)
    % Pris seuls, ils ressemblent à du bruit car ils ne reçoivent qu'une fraction
    % diluée du pool d'origine. Leur force n'apparaît que lorsqu'ils sont sommés.
    for st = 1:num_stimuli
        pooled_rate = num_coll * base_rate * 0.85; % Dilution ajustée pour effondrer la performance unitaire
        approx_spikes = round((t2 - t1) * pooled_rate * 3) + 10;
        uniform_samples = rand(1, approx_spikes);
        intervals = refrac - log(1 - uniform_samples) / pooled_rate;
        spikes_pooled = cumsum(intervals);
        spikes_pooled = spikes_pooled(spikes_pooled >= t1 & spikes_pooled <= t2);
        num_spikes = length(spikes_pooled);
        
        for rp = 1:num_repetitions
            if num_spikes > 0
                shuffled_indices = randperm(num_spikes);
                for nc = 1:num_coll 
                    idx_assigned = shuffled_indices(nc:num_coll:end);
                    CellMatrix{nc, st, rp} = sort(spikes_pooled(idx_assigned));
                end
            else
                for nc = 1:num_coll
                    CellMatrix{nc, st, rp} = [];
                end
            end
        end
    end
    
    %% 2. Generation of INDIVIDUAL coding neurons (Channels: num_coll + 1 to num_coll + num_indi)
    % Chaque neurone possède son propre pattern temporel stable par stimulus.
    % On lui applique un petit "jitter" (gigue temporelle) à chaque essai pour mimer la dégradation de la Fig 7.
    for c_idx = 1:num_indi
        nc = num_coll + c_idx; % Canaux 5 à 7 si num_coll=4
        
        for st = 1:num_stimuli
            % On fige temporairement la graine pour générer le même pattern de base pour ce stimulus
            rng(nc * 150 + st); 
            
            local_rate = base_rate * 1.1; 
            approx_spikes = round((t2 - t1) * local_rate * 3) + 10;
            uniform_samples = rand(1, approx_spikes);
            intervals = refrac - log(1 - uniform_samples) / local_rate;
            spikes_reference = cumsum(intervals);
            spikes_reference = spikes_reference(spikes_reference >= t1 & spikes_reference <= t2);
            
            % On libère la graine pour le reste de la simulation
            rng('shuffle');
            
            for rp = 1:num_repetitions
                if ~isempty(spikes_reference)
                    % Ajout d'une gigue temporelle contrôlée (+/- 8 ms) pour bruiter le signal unitaire
                    jitter = (rand(size(spikes_reference)) - 0.5) * 0.016; 
                    spikes_jittered = sort(spikes_reference + jitter);
                    
                    % Protection des bornes temporelles
                    CellMatrix{nc, st, rp} = spikes_jittered(spikes_jittered >= t1 & spikes_jittered <= t2);
                else
                    CellMatrix{nc, st, rp} = [];
                end
            end
        end
    end
    
    %% 3. Generation of NON-CODING neurons (Channels: num_coll + num_indi + 1 to num_neurons)
    % Bruit de fond uniforme sans aucune structure reproductible ou corrélée.
    for st = 1:num_stimuli
        for rp = 1:num_repetitions
            for nc = (num_coll + num_indi + 1):num_neurons
                approx_spikes = round((t2 - t1) * base_rate * 1.0) + 10;
                uniform_samples = rand(1, approx_spikes);
                intervals = refrac - log(1 - uniform_samples) / base_rate;
                spikes_noise = cumsum(intervals);
                CellMatrix{nc, st, rp} = spikes_noise(spikes_noise >= t1 & spikes_noise <= t2);
            end
        end
    end
    
    %% =========================================================================
    %% SECTION GRAPHIQUE (Raster plots optionnels)
    %% =========================================================================
    if plotting == true
        if other_figs == true
            color_coll = [0.8500 0.3250 0.0980];  % Orange pour Collectif
            color_indi = [1 0 0];                 % Rouge pour Individuel
            color_noise = [0.0000 0.4470 0.7410]; % Bleu pour Bruit (NC)
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
                        if nc <= num_coll
                            current_color = color_coll; line_width = 1.4;
                        elseif nc <= (num_coll + num_indi)
                            current_color = color_indi; line_width = 1.4;
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
            title('Artificial Dataset (Orange: Coll | Red: Indi | Blue: NC Noise)', 'FontSize', 12, 'FontWeight', 'bold');
            hold off;
        end
    end
end