% %% Simulated Annealing algorithm script (for SP hypothesis) 
% % Date: June 2026
% % Author : Laure WOLFF
% 
% function [nb_iterations] = f_simulated_annealing(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, showing, plotting, other_figs)
% %% 1. Initialization of the variables
% cooling_factor = 0.9;        
% alpha_threshold = 1e-5;       
% iterations_per_temp = 5 * num_neurons; 
% N0 = 50;
% max_paliers_est = 200; 
% Matrix_Grid = NaN(max_paliers_est, num_neurons);
% history_perf = zeros(1, max_paliers_est); 
% mask_0 = randi([0, 1], num_neurons, 1);
% if sum(mask_0) == 0, mask_0(randi(num_neurons)) = 1; end
% if sum(mask_0) == num_neurons, mask_0(randi(num_neurons)) = 0; end
% % [P_0,~] = calculate_integrated_P_optimized(CellMatrix, mask_0, ...
% %     num_stimuli, num_repetitions, t1, t2, metric_choice);
% [P_0,~] = calculate_integrated_P_optimized(CellMatrix, mask_0, ...
%     num_stimuli, num_repetitions, t1, t2);
% best_perf_overall = P_0;
% best_mask_overall = mask_0;
% temp_mask = mask_0;
% temp_perf = P_0;
% delta_down = zeros(1, N0);
% count = 0;
% 
% %% 2. Finding T_0
% for n = 1:N0
%     idx = randi(num_neurons);
%     next_mask = temp_mask;
%     next_mask(idx) = 1 - temp_mask(idx);
% 
%     if sum(next_mask) == 0 || sum(next_mask) == num_neurons, continue; end
% 
%     % [next_perf,~] = calculate_integrated_P_optimized(CellMatrix, next_mask, ...
%     %     num_stimuli, num_repetitions, t1, t2, metric_choice);
%      [next_perf,~] = calculate_integrated_P_optimized(CellMatrix, next_mask, ...
%         num_stimuli, num_repetitions, t1, t2);
% 
%     if next_perf <= temp_perf
%         count = count + 1;
%         delta_down(count) = abs(next_perf - temp_perf);
%     end
%     temp_perf = next_perf;
%     temp_mask = next_mask;
% end
% if count > 0
%     % On extrait la partie du tableau remplie
%     filled_deltas = delta_down(1:count);
% 
%     % On ne garde QUE les deltas qui sont des nombres réels finis (on vire -Inf, Inf, NaN)
%     % et non nuls (car un delta de 0 fausse la moyenne thermique)
%     valid_deltas = filled_deltas(isfinite(filled_deltas) & (filled_deltas ~= 0));
% 
%     if ~isempty(valid_deltas)
%         mean_delta = mean(valid_deltas);
%     else
%         mean_delta = 0.005; % Repli (Fallback) si tout le tableau contenait du -Inf
%         if showing
%             fprintf('  [Warning SA] Tous les tirages initiaux ont touché le neurone piège (-Inf). Default delta appliqué.\n');
%         end
%     end
% else
%     mean_delta = 0.005;
% end
% 
% % Calcul de T_0
% T_0 = - mean_delta / log(0.95); 
% 
% % Ultime barrière de sécurité pour empêcher le crash de la fonction
% if T_0 <= 1e-7 || isnan(T_0) || isinf(T_0)
%     T_0 = 0.5; % On force une température de démarrage standard (0.5) pour sauver le script
%     if showing
%         fprintf('  [Warning SA] T_0 calculé invalide. Température forcée à %.2f pour éviter le crash.\n', T_0);
%     end
% end
% if showing, fprintf('T_0 found: %.6f \n', T_0); end
% max_iter_est = max_paliers_est * iterations_per_temp;
% hist_iter_P      = zeros(1, max_iter_est);
% hist_iter_bestP  = zeros(1, max_iter_est);
% hist_iter_size   = zeros(1, max_iter_est);
% hist_iter_temp   = zeros(1, max_iter_est);
% 
% %% 3. Simulated Annealing Loop
% theta = T_0;            
% unchanged_temp_cycles = 0;
% palier_idx = 0;
% nb_iterations = 0;
% while theta > alpha_threshold
%     palier_idx = palier_idx + 1;
% 
%     % Security allocation
%     if palier_idx > size(Matrix_Grid, 1)
%         % If limit is exceeded (highly rare), double the size at once to minimize copies
%         new_size = size(Matrix_Grid, 1) * 2;
% 
%         Matrix_Grid_Expanded = NaN(new_size, num_neurons);
%         Matrix_Grid_Expanded(1:size(Matrix_Grid, 1), :) = Matrix_Grid;
%         Matrix_Grid = Matrix_Grid_Expanded;
% 
%         history_perf_Expanded = zeros(1, new_size);
%         history_perf_Expanded(1:length(history_perf)) = history_perf;
%         history_perf = history_perf_Expanded;
%     end
% 
%     if showing 
%         fprintf('Temp: %.6f | Current P: %.4f\n', theta, temp_perf); 
%     end
% 
%     for iter = 1:iterations_per_temp
%         nb_iterations= nb_iterations+1;
%         active_count = sum(temp_mask);
%         next_mask = temp_mask;
% 
%         % Security
%         if active_count == 1
%             zero_indices = find(temp_mask == 0);
%             idx_explore = zero_indices(randi(length(zero_indices)));
%             next_mask(idx_explore) = 1;
%         elseif active_count == num_neurons
%             one_indices = find(temp_mask == 1);
%             idx_explore = one_indices(randi(length(one_indices)));
%             next_mask(idx_explore) = 0;
%         else
%             idx_explore = randi(num_neurons);
%             next_mask(idx_explore) = 1 - temp_mask(idx_explore);
%         end
% 
%         % [next_perf,~] = calculate_integrated_P_optimized(CellMatrix, next_mask, ...
%         %     num_stimuli, num_repetitions, t1, t2, metric_choice);
%         [next_perf,~] = calculate_integrated_P_optimized(CellMatrix, next_mask, ...
%             num_stimuli, num_repetitions, t1, t2);
% 
%         if next_perf > temp_perf
%             temp_mask = next_mask;
%             temp_perf = next_perf;
%         else
%             q = exp(-abs(next_perf - temp_perf) / theta); 
%             if rand() < q
%                 temp_mask = next_mask;
%                 temp_perf = next_perf;
%             end
%         end
% 
%         if temp_perf > best_perf_overall
%             best_perf_overall = temp_perf;
%             best_mask_overall = temp_mask;
%         end
%         hist_iter_P(nb_iterations)     = temp_perf;
%         hist_iter_bestP(nb_iterations) = best_perf_overall;
%         hist_iter_size(nb_iterations)  = sum(temp_mask);
%         hist_iter_temp(nb_iterations)  = theta;
%     end
%     Matrix_Grid(palier_idx, :) = temp_mask';
%     history_perf(palier_idx) = temp_perf;
% 
%     if palier_idx >= 2 && abs(history_perf(palier_idx) - history_perf(palier_idx-1)) < 1e-6
%         % unchanged_temp_cycles = unchanged_temp_cycles + 1;
%         % if unchanged_temp_cycles == 2
%             if showing
%                 fprintf(['Exit: Performance remained unchanged for 2 ' ...
%                     'consecutive temperature cycles.\n']);
%             end
%             break; 
%     end
% 
%     theta = theta * cooling_factor;
% end
% Matrix_Grid = Matrix_Grid(1:palier_idx, :);
% history_perf = history_perf(1:palier_idx);
% hist_iter_P     = hist_iter_P(1:nb_iterations);
% hist_iter_bestP = hist_iter_bestP(1:nb_iterations);
% hist_iter_size  = hist_iter_size(1:nb_iterations);
% hist_iter_temp  = hist_iter_temp(1:nb_iterations);
% 
% %% 4. Final Wrap-up
% best_subpop = find(best_mask_overall == 1)';
% if showing
%     fprintf('Optimal subpopulation found: [%s]\n', num2str(best_subpop));
%     fprintf('Max performance P = %.4f\n', best_perf_overall);
%     fprintf('number of iteration %.4f\n', nb_iterations);
% end
% 
% %% 5. Plotting (With Dynamic Adaptive Scales)
% if plotting == true && ~isempty(Matrix_Grid)
%     num_paliers_reals = size(Matrix_Grid, 1);
% 
%     if num_neurons <= 15
%         tick_step_X = 1;
%     elseif num_neurons <= 40
%         tick_step_X = 5;
%     else
%         tick_step_X = 10; % Un tick tous les 10 neurones si N=100
%     end
% 
%     if num_paliers_reals <= 20
%         tick_step_Y = 1;
%     elseif num_paliers_reals <= 60
%         tick_step_Y = 5;
%     else
%         tick_step_Y = 10; % Un tick tous les 10 paliers si la simulation est longue
%     end
% 
%     figure('Name', 'Results - Simulated Annealing','Color', [1 1 1]);
% 
%     % Matrix which shows the several masks
%     subplot(1, 4, 1:3);
%     imagesc(1:num_neurons, 1:num_paliers_reals, Matrix_Grid);
%     mymap = [0.2 0.4 0.8; 0.9 0.2 0.2]; 
%     colormap(gca, mymap); 
%     clim([0, 1]);
%     set(gca, 'YDir', 'normal'); 
%     hold on;
% 
%     % Adding cross to mark the best subpopulation (at the top row)
%     for i = 1:length(best_subpop)
%         plot(best_subpop(i), num_paliers_reals, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
%     end
% 
%     box on;
%     set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 10);
%     set(gca, 'XTick', 1:tick_step_X:num_neurons);
%     set(gca, 'YTick', 1:tick_step_Y:num_paliers_reals);
% 
%     xlabel('Neuron ID', 'FontSize', 12, 'FontWeight', 'bold');
%     ylabel('Temperature Steps (Cooling)', 'FontSize', 12, 'FontWeight', 'bold');
%     title('Simulated Annealing - Selected Neurons History', 'FontSize', 13, 'FontWeight', 'bold');
% 
%     % Plot's legends
%     cb = colorbar;
%     set(cb, 'Ticks', [0.25, 0.75], 'TickLabels', {'Desactivated (0)', 'Activated (1)'}, 'FontSize', 10);
% 
%     % Performance tracking plot (Right subplot)
%     subplot(1, 4, 4);
%     plot(history_perf, 1:num_paliers_reals, '-ko', 'LineWidth', 1.5, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 4);
%     hold on;
% 
%     plot(best_perf_overall, num_paliers_reals, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [1 1 1]);
%     plot(best_perf_overall, num_paliers_reals, 'rx', 'MarkerSize', 6, 'LineWidth', 1.5);
% 
%     grid on; box on;
%     set(gca, 'YDir', 'normal', 'TickDir', 'out', 'LineWidth', 1.2, 'YTickLabel', [], 'FontSize', 10);
% 
%     set(gca, 'YTick', 1:tick_step_Y:num_paliers_reals);
%     ylim([0.5, num_paliers_reals + 0.5]);
%     xlabel('Performance P', 'FontSize', 12, 'FontWeight', 'bold');
%     title('P(temp)', 'FontSize', 13, 'FontWeight', 'bold');
% 
%     hold off; 
% 
%     %% Convergence Diagnostics Figure
%     if other_figs == true
%         figure('Name', 'Simulated Annealing - Convergence Diagnostics', 'Color', 'w', 'Position', [150, 150, 1000, 400]);
%         tiledlayout(1, 3, 'TileSpacing', 'compact');
% 
%         % 1. Discrimination performance
%         nexttile
%         plot(hist_iter_P, 'k', 'LineWidth', 1.2);
%         hold on;
%         plot(hist_iter_bestP, 'r', 'LineWidth', 2);
%         xlabel('Iteration');
%         ylabel('Performance P');
%         legend({'Current', 'Best'}, 'Location', 'best');
%         title('Discrimination performance');
%         grid on; box on;
% 
%         % 2. Population size
%         nexttile
%         plot(hist_iter_size, 'b', 'LineWidth', 1.2);
%         hold on;
%         yline(length(best_subpop), '--r', 'Optimal found', 'LineWidth', 1.5);
%         xlabel('Iteration');
%         ylabel('Population size');
%         title('Subpopulation Size Track');
%         grid on; box on;
% 
%         % 3. Temperature cooling 
%         nexttile
%         semilogy(hist_iter_temp, 'm', 'LineWidth', 1.5);
%         xlabel('Iteration');
%         ylabel('Temperature');
%         title('Cooling Schedule (Log Scale)');
%         grid on; box on;
% 
%         % Global title
%         sgtitle(sprintf('Simulated Annealing Dynamics | Best P = %.4f | Best Population = [%s]',...
%             best_perf_overall, num2str(sort(best_subpop))));
%     end
%     shg;
% end
% end


%% Simulated Annealing algorithm script (for SP hypothesis with memoization cache) 
% Date: July 2026
% Author : Laure WOLFF

function [nb_iterations] = f_simulated_annealing(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, showing, plotting, other_figs)
%% 0. Memoization Cache Initialization
% containers.Map holds (Key: String Mask, Value: Performance P)
memo_cache = containers.Map('KeyType', 'char', 'ValueType', 'double');
cache_hits = 0; % Counter to track saved evaluations

%% 1. Initialization of the variables
cooling_factor = 0.9;        
alpha_threshold = 1e-5;       
iterations_per_temp = 5 * num_neurons; 
N0 = 50;
max_paliers_est = 200; 
Matrix_Grid = NaN(max_paliers_est, num_neurons);
history_perf = zeros(1, max_paliers_est); 

mask_0 = randi([0, 1], num_neurons, 1);
if sum(mask_0) == 0, mask_0(randi(num_neurons)) = 1; end
if sum(mask_0) == num_neurons, mask_0(randi(num_neurons)) = 0; end

% First evaluation using memoized function
[P_0, memo_cache, cache_hits] = evaluate_mask_memoized(...
    CellMatrix, mask_0, num_stimuli, num_repetitions, t1, t2, memo_cache, cache_hits);

best_perf_overall = P_0;
best_mask_overall = mask_0;
temp_mask = mask_0;
temp_perf = P_0;
delta_down = zeros(1, N0);
count = 0;

%% 2. Finding T_0
for n = 1:N0
    idx = randi(num_neurons);
    next_mask = temp_mask;
    next_mask(idx) = 1 - temp_mask(idx);
    
    if sum(next_mask) == 0 || sum(next_mask) == num_neurons, continue; end
    
    [next_perf, memo_cache, cache_hits] = evaluate_mask_memoized(...
        CellMatrix, next_mask, num_stimuli, num_repetitions, t1, t2, memo_cache, cache_hits);
        
    if next_perf <= temp_perf
        count = count + 1;
        delta_down(count) = abs(next_perf - temp_perf);
    end
    temp_perf = next_perf;
    temp_mask = next_mask;
end

if count > 0
    filled_deltas = delta_down(1:count);
    valid_deltas = filled_deltas(isfinite(filled_deltas) & (filled_deltas ~= 0));
    
    if ~isempty(valid_deltas)
        mean_delta = mean(valid_deltas);
    else
        mean_delta = 0.005; % Fallback
        if showing
            fprintf('  [Warning SA] All initial draws hit trap neurons (-Inf). Default delta applied.\n');
        end
    end
else
    mean_delta = 0.005;
end

% Temperature T_0 calculation
T_0 = - mean_delta / log(0.95); 
if T_0 <= 1e-7 || isnan(T_0) || isinf(T_0)
    T_0 = 0.5; % Safety guard
    if showing
        fprintf('  [Warning SA] Calculated T_0 invalid. Forced to %.2f.\n', T_0);
    end
end
if showing, fprintf('T_0 found: %.6f \n', T_0); end

max_iter_est = max_paliers_est * iterations_per_temp;
hist_iter_P      = zeros(1, max_iter_est);
hist_iter_bestP  = zeros(1, max_iter_est);
hist_iter_size   = zeros(1, max_iter_est);
hist_iter_temp   = zeros(1, max_iter_est);

%% 3. Simulated Annealing Loop
theta = T_0;            
palier_idx = 0;
nb_iterations = 0;

while theta > alpha_threshold
    palier_idx = palier_idx + 1;
    
    % Dynamic expansion if limit exceeded
    if palier_idx > size(Matrix_Grid, 1)
        new_size = size(Matrix_Grid, 1) * 2;
        
        Matrix_Grid_Expanded = NaN(new_size, num_neurons);
        Matrix_Grid_Expanded(1:size(Matrix_Grid, 1), :) = Matrix_Grid;
        Matrix_Grid = Matrix_Grid_Expanded;
        
        history_perf_Expanded = zeros(1, new_size);
        history_perf_Expanded(1:length(history_perf)) = history_perf;
        history_perf = history_perf_Expanded;
    end
    
    if showing 
        fprintf('Temp: %.6f | Current P: %.4f | Cache Hits: %d\n', theta, temp_perf, cache_hits); 
    end
    
    for iter = 1:iterations_per_temp
        nb_iterations = nb_iterations + 1;
        active_count = sum(temp_mask);
        next_mask = temp_mask;
        
        % Neighborhood exploration
        if active_count == 1
            zero_indices = find(temp_mask == 0);
            idx_explore = zero_indices(randi(length(zero_indices)));
            next_mask(idx_explore) = 1;
        elseif active_count == num_neurons
            one_indices = find(temp_mask == 1);
            idx_explore = one_indices(randi(length(one_indices)));
            next_mask(idx_explore) = 0;
        else
            idx_explore = randi(num_neurons);
            next_mask(idx_explore) = 1 - temp_mask(idx_explore);
        end
        
        % Memoized Performance Evaluation
        [next_perf, memo_cache, cache_hits] = evaluate_mask_memoized(...
            CellMatrix, next_mask, num_stimuli, num_repetitions, t1, t2, memo_cache, cache_hits);
        
        if next_perf > temp_perf
            temp_mask = next_mask;
            temp_perf = next_perf;
        else
            q = exp(-abs(next_perf - temp_perf) / theta); 
            if rand() < q
                temp_mask = next_mask;
                temp_perf = next_perf;
            end
        end
        
        if temp_perf > best_perf_overall
            best_perf_overall = temp_perf;
            best_mask_overall = temp_mask;
        end
        hist_iter_P(nb_iterations)     = temp_perf;
        hist_iter_bestP(nb_iterations) = best_perf_overall;
        hist_iter_size(nb_iterations)  = sum(temp_mask);
        hist_iter_temp(nb_iterations)  = theta;
    end
    
    Matrix_Grid(palier_idx, :) = temp_mask';
    history_perf(palier_idx) = temp_perf;
    
    if palier_idx >= 2 && abs(history_perf(palier_idx) - history_perf(palier_idx-1)) < 1e-6
        if showing
            fprintf('Exit: Performance remained unchanged for 2 consecutive temperature cycles.\n');
        end
        break; 
    end
    
    theta = theta * cooling_factor;
end

Matrix_Grid = Matrix_Grid(1:palier_idx, :);
history_perf = history_perf(1:palier_idx);
hist_iter_P     = hist_iter_P(1:nb_iterations);
hist_iter_bestP = hist_iter_bestP(1:nb_iterations);
hist_iter_size  = hist_iter_size(1:nb_iterations);
hist_iter_temp  = hist_iter_temp(1:nb_iterations);

%% 4. Final Wrap-up
best_subpop = find(best_mask_overall == 1)';
if showing
    fprintf('\n--- SIMULATED ANNEALING SUMMARY ---\n');
    fprintf('Optimal subpopulation found: [%s]\n', num2str(best_subpop));
    fprintf('Max performance P = %.4f\n', best_perf_overall);
    fprintf('Total iterations: %d\n', nb_iterations);
    fprintf('Unique masks evaluated: %d\n', memo_cache.Count);
    fprintf('Cache hits (Calculations saved!): %d (%.1f%% saved)\n', ...
        cache_hits, (cache_hits / (nb_iterations + N0 + 1)) * 100);
end

%% 5. Plotting (With Dynamic Adaptive Scales)
if plotting == true && ~isempty(Matrix_Grid)
    num_paliers_reals = size(Matrix_Grid, 1);
    
    if num_neurons <= 15
        tick_step_X = 1;
    elseif num_neurons <= 40
        tick_step_X = 5;
    else
        tick_step_X = 10;
    end
    
    if num_paliers_reals <= 20
        tick_step_Y = 1;
    elseif num_paliers_reals <= 60
        tick_step_Y = 5;
    else
        tick_step_Y = 10;
    end
    
    figure('Name', 'Results - Simulated Annealing', 'Color', [1 1 1]);
    
    subplot(1, 4, 1:3);
    imagesc(1:num_neurons, 1:num_paliers_reals, Matrix_Grid);
    mymap = [0.2 0.4 0.8; 0.9 0.2 0.2]; 
    colormap(gca, mymap); 
    clim([0, 1]);
    set(gca, 'YDir', 'normal'); 
    hold on;
    
    for i = 1:length(best_subpop)
        plot(best_subpop(i), num_paliers_reals, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
    
    box on;
    set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 10);
    set(gca, 'XTick', 1:tick_step_X:num_neurons);
    set(gca, 'YTick', 1:tick_step_Y:num_paliers_reals);
    
    xlabel('Neuron ID', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Temperature Steps (Cooling)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Simulated Annealing - Selected Neurons History', 'FontSize', 13, 'FontWeight', 'bold');
    
    cb = colorbar;
    set(cb, 'Ticks', [0.25, 0.75], 'TickLabels', {'Desactivated (0)', 'Activated (1)'}, 'FontSize', 10);
    
    subplot(1, 4, 4);
    plot(history_perf, 1:num_paliers_reals, '-ko', 'LineWidth', 1.5, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 4);
    hold on;
    
    plot(best_perf_overall, num_paliers_reals, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [1 1 1]);
    plot(best_perf_overall, num_paliers_reals, 'rx', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    grid on; box on;
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'LineWidth', 1.2, 'YTickLabel', [], 'FontSize', 10);
    
    set(gca, 'YTick', 1:tick_step_Y:num_paliers_reals);
    ylim([0.5, num_paliers_reals + 0.5]);
    xlabel('Performance P', 'FontSize', 12, 'FontWeight', 'bold');
    title('P(temp)', 'FontSize', 13, 'FontWeight', 'bold');
    
    hold off; 
    
    if other_figs == true
        figure('Name', 'Simulated Annealing - Convergence Diagnostics', 'Color', 'w', 'Position', [150, 150, 1000, 400]);
        tiledlayout(1, 3, 'TileSpacing', 'compact');
        
        nexttile
        plot(hist_iter_P, 'k', 'LineWidth', 1.2);
        hold on;
        plot(hist_iter_bestP, 'r', 'LineWidth', 2);
        xlabel('Iteration');
        ylabel('Performance P');
        legend({'Current', 'Best'}, 'Location', 'best');
        title('Discrimination performance');
        grid on; box on;
        
        nexttile
        plot(hist_iter_size, 'b', 'LineWidth', 1.2);
        hold on;
        yline(length(best_subpop), '--r', 'Optimal found', 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('Population size');
        title('Subpopulation Size Track');
        grid on; box on;
        
        nexttile
        semilogy(hist_iter_temp, 'm', 'LineWidth', 1.5);
        xlabel('Iteration');
        ylabel('Temperature');
        title('Cooling Schedule (Log Scale)');
        grid on; box on;
        
        sgtitle(sprintf('Simulated Annealing Dynamics | Best P = %.4f | Best Population = [%s]',...
            best_perf_overall, num2str(sort(best_subpop))));
    end
    shg;
end
end


%% --- HELPER FUNCTION FOR MEMOIZED EVALUATION ---
function [perf, memo_cache, cache_hits] = evaluate_mask_memoized(CellMatrix, mask, S, R, t1, t2, memo_cache, cache_hits)
    % Convert mask array to unique string key (e.g., [1; 0; 1] -> '101')
    key = sprintf('%d', mask(:));
    
    if isKey(memo_cache, key)
        % Retrieve cached value (Instantaneous 0 ms!)
        perf = memo_cache(key);
        cache_hits = cache_hits + 1;
    else
        % Compute full performance function
        [perf, ~] = calculate_integrated_P_optimized(CellMatrix, mask, S, R, t1, t2);
        % Save to cache for future encounters
        memo_cache(key) = perf;
    end
end