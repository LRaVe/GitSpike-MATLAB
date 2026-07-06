%% Script to calculate and plot the performance matrices (Figures 4E, F, G)
% author: Laure WOLFF
% Date: July 2026
function [P_all_neurons, P_pop, M_max, opt_LL, PLL_total] = calculate_and_plot_performance_matrix(All_Matrices_M, All_Matrix_D, num_neurons, num_stimuli, num_repetitions, plotting)
    P_all_neurons = zeros(num_stimuli, num_stimuli, num_neurons);
    
    % --- BLOCK 1: INDIVIDUAL NEURON PERFORMANCE (FIGURE 4E) ---
    for n = 1:num_neurons
        MatrixM = All_Matrices_M(:, :, n); 
        distance_matrix = All_Matrix_D(:, :, n);
        P_matrix_current = zeros(num_stimuli, num_stimuli);
        
        % Extract maximum distance to normalize each neuron profile dynamically
        max_dist_neuron = max(distance_matrix(:));
        if max_dist_neuron == 0, max_dist_neuron = 1; end
        
        for st1 = 1:num_stimuli
            idx1 = (st1-1)*num_repetitions + 1 : st1*num_repetitions;
            for st2 = 1:num_stimuli
                if st1 == st2, continue; end
                
                % If Wilcoxon test confirmed statistical discrimination
                if MatrixM(st1, st2) == 1 
                    idx2 = (st2-1)*num_repetitions + 1 : st2*num_repetitions;
                    mean_inter = mean(reshape(distance_matrix(idx1, idx2), [], 1));
                    
                    % Compute upper triangular intra-distances to avoid self-comparison bias
                    intra_1 = distance_matrix(idx1, idx1); dist_intra_1 = mean(intra_1(triu(true(num_repetitions), 1)));
                    intra_2 = distance_matrix(idx2, idx2); dist_intra_2 = mean(intra_2(triu(true(num_repetitions), 1)));
                    
                    contrast = mean_inter - (dist_intra_1 + dist_intra_2)/2;
                    
                    % Safe dynamic rescaling bounded between 0 and 1
                    if contrast > 0
                        P_matrix_current(st1, st2) = min(1, (contrast / max_dist_neuron) * 2.5);
                    else
                        P_matrix_current(st1, st2) = 0;
                    end
                else
                    % Failed Wilcoxon test -> Absolute zero (Deep Blue)
                    P_matrix_current(st1, st2) = 0;
                end
            end
        end
        P_all_neurons(:, :, n) = P_matrix_current;
    end
    
    % --- BLOCK 2: POPULATION PERFORMANCE MATRIX  ---
    P_pop = max(P_all_neurons, [], 3);
    
    % --- BLOCK 3: BEST NEURON MAPPING  ---
    M_max = zeros(num_stimuli, num_stimuli);
    
    % Global selectivity count (total active neurons in Wilcoxon mask)
    global_selectivity = zeros(num_neurons, 1);
    for n = 1:num_neurons
        global_selectivity(n) = sum(sum(All_Matrices_M(:, :, n)));
    end

    for s = 1:num_stimuli
        for s_prime = 1:num_stimuli
            if s == s_prime, continue; end
            
            if P_pop(s, s_prime) > 0
                perf_profile = squeeze(P_all_neurons(s, s_prime, :));
                is_valid_coder = squeeze(All_Matrices_M(s, s_prime, :));
                
                % Z-score normalization of profiles to equalize firing rate differences
                mean_p = mean(perf_profile(is_valid_coder == 1));
                std_p = std(perf_profile(is_valid_coder == 1));
                if isempty(std_p) || std_p == 0, std_p = 1; end
                
                normalized_profile = (perf_profile - mean_p) ./ std_p;
                
                % Find candidates among valid Wilcoxon coders
                max_norm_perf = max(normalized_profile .* is_valid_coder);
                candidates = find(normalized_profile >= (max_norm_perf - 0.1) & is_valid_coder == 1);
                
                if ~isempty(candidates)
                    % Priority rule: pick the candidate with the fewest active cells overall
                    candidate_scores = global_selectivity(candidates);
                    [~, rel_idx] = min(candidate_scores);
                    M_max(s, s_prime) = candidates(rel_idx);
                else
                    % Standard fallback
                    [~, backup_idx] = max(perf_profile .* is_valid_coder);
                    if backup_idx > 0 && is_valid_coder(backup_idx) == 1
                        M_max(s, s_prime) = backup_idx;
                    else
                        M_max(s, s_prime) = 0;
                    end
                end
            else
                M_max(s, s_prime) = 0; 
            end
        end
    end
    
    % --- BLOCK 4: OPTIMIZED POPULATION SELECTION (EQ 14) ---
    opt_LL = unique(M_max(M_max > 0))';
    
    % --- BLOCK 5: GLOBAL PERFORMANCE SCORES (EQ 15) ---
    mask_distinct_pairs = ~eye(num_stimuli);
    PLL_total = mean(P_pop(mask_distinct_pairs));
    
    % --- PLOTTING SECTION ---
    if plotting == true
        % Figure 4E
        figure('Name', 'Individual performance matrices P_n', 'Color', 'w');
        cols = ceil(sqrt(num_neurons)); rows = ceil(num_neurons / cols);
        for n = 1:num_neurons
            subplot(rows, cols, n);
            imagesc(P_all_neurons(:, :, n));
            colormap(gca, jet); clim([0 1]); colorbar; axis square;
            title(sprintf('P_{%d} (Neuron %d)', n, n), 'FontSize', 10, 'FontWeight', 'bold');
            set(gca, 'XTick', 1:num_stimuli, 'YTick', 1:num_stimuli);
        end
        
        % Figures 4F & 4G
        figure('Name', 'Figures 4F & 4G: Population discrimination analysis', 'Color', 'w');
        
        % Fig 4F
        subplot(1, 2, 1);
        imagesc(P_pop); colormap(gca, jet); clim([0 1]); colorbar; axis square;
        set(gca, 'XTick', 1:num_stimuli, 'YTick', 1:num_stimuli);
        title('P_{max} (Population Performance)', 'FontWeight', 'bold');
        subtitle(sprintf('Global LL performance | P_{LL} = %.4f', PLL_total));
        
        % --- Fig 4G Plot ---
        subplot(1, 2, 2);
        imagesc(M_max); 
        
        % Base palette for the first 4 neurons
        custom_colors = [
            0.0, 0.0, 0.0;   % 0 -> Pure Black (No active winner)
            0.0, 0.45, 1.0;  % 1 -> Pure Blue (Neuron 1)
            1.0, 0.0, 0.0;   % 2 -> Pure Red (Neuron 2)
            0.0, 0.65, 0.0;  % 3 -> Pure Green (Neuron 3)
            1.0, 0.85, 0.0;  % 4 -> Bright Yellow (Neuron 4)
        ];
        
        if (num_neurons + 1) > size(custom_colors, 1)
            num_extra = (num_neurons + 1) - size(custom_colors, 1);
            % Using lines() or jet() to generate high-contrast extra colors
            extra_map = lines(num_extra); 
            custom_colors = [custom_colors; extra_map];
        end
        
        % Apply the perfectly sized colormap
        colormap(gca, custom_colors(1:num_neurons + 1, :));
        clim([0 num_neurons]); 
        
        cb = colorbar; 
        set(cb, 'Ticks', 0:num_neurons, 'TickLabels', 0:num_neurons);
        axis square;
        set(gca, 'XTick', 1:num_stimuli, 'YTick', 1:num_stimuli);
        title('M_{max} (Best Neuron Mapping - Fig 4G)', 'FontWeight', 'bold');
        subtitle(sprintf('Optimized Subpopulation opt_{LL} = [%s]', num2str(opt_LL)));
    end
end