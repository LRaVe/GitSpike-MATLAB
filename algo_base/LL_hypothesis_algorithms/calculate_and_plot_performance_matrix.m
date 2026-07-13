%% Script to calculate and plot the performance matrices (Figures 4E, F, G)
% Author: Laure WOLFF
% Date: July 2026
function [P_all_neurons, P_pop, M_max, opt_LL, PLL_total] = calculate_and_plot_performance_matrix(All_Matrices_M, All_Matrix_D, num_neurons, num_stimuli, num_repetitions, plotting)

    P_all_neurons = zeros(num_stimuli, num_stimuli, num_neurons);
    mask_distinct_pairs = ~eye(num_stimuli);

% --- BLOCK 1: INDIVIDUAL NEURON PERFORMANCE---
    for n = 1:num_neurons
        MatrixM = All_Matrices_M(:, :, n); 
        distance_matrix = All_Matrix_D(:, :, n);
        P_matrix_current = zeros(num_stimuli, num_stimuli);

        for st1 = 1:num_stimuli
            idx1 = (st1-1)*num_repetitions + 1 : st1*num_repetitions;
            for st2 = 1:num_stimuli
                if st1 == st2, continue; end

                if MatrixM(st1, st2) == 0
                    P_matrix_current(st1, st2) = 0;
                    continue;
                end

                idx2 = (st2-1)*num_repetitions + 1 : st2*num_repetitions;

                mean_inter = mean(reshape(distance_matrix(idx1, idx2), [], 1));

                intra_1 = distance_matrix(idx1, idx1); 
                dist_intra_1 = mean(intra_1(triu(true(num_repetitions), 1)));

                intra_2 = distance_matrix(idx2, idx2); 
                dist_intra_2 = mean(intra_2(triu(true(num_repetitions), 1)));

                contrast = mean_inter - (dist_intra_1 + dist_intra_2)/2;

                if contrast > 0
                    P_matrix_current(st1, st2) = contrast / (mean_inter + dist_intra_1 + dist_intra_2);
                else
                    P_matrix_current(st1, st2) = 0;
                end
            end
        end
        P_all_neurons(:, :, n) = P_matrix_current;
    end

    % --- BLOCK 2: POPULATION PERFORMANCE MATRIX ---
    P_pop = max(P_all_neurons, [], 3);
    PLL_total = mean(P_pop(mask_distinct_pairs));

    % --- BLOCK 3: BEST NEURON MAPPING & TIE-BREAKER ---
    M_max = zeros(num_stimuli, num_stimuli);

    global_selectivity = squeeze(sum(sum(All_Matrices_M, 1), 2));

    for s = 1:num_stimuli
        for s_prime = 1:num_stimuli
            if s == s_prime || P_pop(s, s_prime) == 0, continue; end

            perf_profile = squeeze(P_all_neurons(s, s_prime, :));
            is_valid_coder = squeeze(All_Matrices_M(s, s_prime, :));

            max_perf = max(perf_profile .* is_valid_coder);
            if max_perf == 0, continue; end

            candidates = find(perf_profile >= (max_perf * 0.95) & is_valid_coder == 1);

            if ~isempty(candidates)
                % Tie-breaker : prendre le neurone le plus sélectif globalement
                [~, rel_idx] = min(global_selectivity(candidates));
                M_max(s, s_prime) = candidates(rel_idx);
            end
        end
    end

    % --- BLOCK 4: OPTIMIZED POPULATION SELECTION ---
    opt_LL = unique(M_max(M_max > 0))';

    % --- PLOTTING SECTION ---
    if plotting == true
        % Figure 4E : Performances individuelles
        figure('Name', 'Individual performance matrices P_n', 'Color', 'w');
        cols = ceil(sqrt(num_neurons)); rows = ceil(num_neurons / cols);
        for n = 1:num_neurons
            subplot(rows, cols, n); 
            imagesc(P_all_neurons(:, :, n));
            colormap(gca, jet); 
            colorbar; 
            axis square;
            title(sprintf('P_{%d} (Neuron %d)', n, n), 'FontSize', 10, 'FontWeight', 'bold');
            set(gca, 'XTick', 1:num_stimuli, 'YTick', 1:num_stimuli);
        end

        % Figures 4F & 4G : Analyse de Population
        figure('Name', 'Population discrimination analysis', 'Color', 'w');

        % Subplot F : P_max
        subplot(1, 2, 1); 
        imagesc(P_pop); 
        colormap(gca, jet);  
        colorbar; 
        axis square;
        set(gca, 'XTick', 1:num_stimuli, 'YTick', 1:num_stimuli);
        title('P_{max} (Population Performance)', 'FontWeight', 'bold');
        subtitle(sprintf('Global LL performance | P_{LL} = %.4f', PLL_total));

        % Subplot G : M_max
        subplot(1, 2, 2); 
        imagesc(M_max); 
        custom_colors = [0,0,0; 0,0.45,1; 1,0,0; 0,0.65,0; 1,0.85,0]; 
        if (num_neurons + 1) > size(custom_colors, 1)
            custom_colors = [custom_colors; lines((num_neurons + 1) - size(custom_colors, 1))];
        end
        colormap(gca, custom_colors(1:num_neurons + 1, :)); 
        clim([0 num_neurons]); 
        cb = colorbar; 
        set(cb, 'Ticks', 0:num_neurons, 'TickLabels', 0:num_neurons); 
        axis square;
        set(gca, 'XTick', 1:num_stimuli, 'YTick', 1:num_stimuli);
        title('M_{max} (Best Neuron Mapping)', 'FontWeight', 'bold');
        subtitle(sprintf('Optimized Subpopulation opt_{LL} = [%s]', num2str(opt_LL)));
    end
end
