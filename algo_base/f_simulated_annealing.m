%% Simulated Annealing algorithm script (for SP hypothesis) 
% Date: June 2026
% Author : Laure WOLFF

function [nb_iterations] = f_simulated_annealing(CellMatrix, num_neurons, num_stimuli, num_repetitions, t1, t2, metric_choice, showing, plotting, other_figs)
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

[P_0,~] = calculate_integrated_P_optimized(CellMatrix, mask_0, ...
    num_stimuli, num_repetitions, t1, t2, metric_choice);
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
    
    [next_perf,~] = calculate_integrated_P_optimized(CellMatrix, next_mask, ...
        num_stimuli, num_repetitions, t1, t2, metric_choice);
        
    if next_perf <= temp_perf
        count = count + 1;
        delta_down(count) = abs(next_perf - temp_perf);
    end
    temp_perf = next_perf;
    temp_mask = next_mask;
end

if count > 0
    mean_delta = mean(delta_down(1:count));
else
    mean_delta = 0.005;
end
T_0 = - mean_delta / log(0.95); 

if T_0 <= 1e-7 || isnan(T_0)
    error("f_simulated_annealing:InvalidTemperature", ...
          "The initial temperature T_0 is null, too small or NaN.");
end
if showing, fprintf('T_0 found: %.6f \n', T_0); end

max_iter_est = max_paliers_est * iterations_per_temp;
hist_iter_P      = zeros(1, max_iter_est);
hist_iter_bestP  = zeros(1, max_iter_est);
hist_iter_size   = zeros(1, max_iter_est);
hist_iter_temp   = zeros(1, max_iter_est);

%% 3. Simulated Annealing Loop
theta = T_0;            
unchanged_temp_cycles = 0;
palier_idx = 0;
nb_iterations = 0;

while theta > alpha_threshold
    palier_idx = palier_idx + 1;
    
    % Security allocation
    if palier_idx > size(Matrix_Grid, 1)
        Matrix_Grid = [Matrix_Grid; NaN(50, num_neurons)]; %#ok<AGROW>
        history_perf = [history_perf, zeros(1, 50)]; %#ok<AGROW>
    end
    
    if showing 
        fprintf('Temp: %.6f | Current P: %.4f\n', theta, temp_perf); 
    end
    
    for iter = 1:iterations_per_temp
        nb_iterations= nb_iterations+1;
        active_count = sum(temp_mask);
        next_mask = temp_mask;
        
        % Security
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
        
        [next_perf,~] = calculate_integrated_P_optimized(CellMatrix, next_mask, ...
            num_stimuli, num_repetitions, t1, t2, metric_choice);
        
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
        unchanged_temp_cycles = unchanged_temp_cycles + 1;
        if unchanged_temp_cycles >= 2
            if showing
                fprintf(['Exit: Performance remained unchanged for 2 ' ...
                    'consecutive temperature cycles.\n']);
            end
            break; 
        end
    else
        unchanged_temp_cycles = 0;
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
    fprintf('Optimal subpopulation found: [%s]\n', num2str(best_subpop));
    fprintf('Max performance P = %.4f\n', best_perf_overall);
    fprintf('number of iteration %.4f\n', nb_iterations);
end

%% 5. Plotting 
if plotting == true && ~isempty(Matrix_Grid)
    num_paliers_reals = size(Matrix_Grid, 1);
    
    figure('Name', 'Results - Simulated Annealing','Color', [1 1 1]);
    
    % Matrix which shox the several mask
    subplot(1, 4, 1:3);
    imagesc(1:num_neurons, 1:num_paliers_reals, Matrix_Grid);
    mymap = [0.2 0.4 0.8; 0.9 0.2 0.2]; 
    colormap(gca, mymap); 
    clim([0, 1]);
    set(gca, 'YTick', 1:num_paliers_reals)
    set(gca, 'YDir', 'normal'); 
    hold on;
    
    % Adding cross to mark the best subpopulation
    for i = 1:length(best_subpop)
        plot(best_subpop(i), num_paliers_reals, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    end
    
    box on;
    set(gca, 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 10);
    set(gca, 'XTick', 1:num_neurons);
    xlabel('Neuron ID', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Temperature Steps (Cooling)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Simulated Annealing - Selected Neurons History', 'FontSize', 13, 'FontWeight', 'bold');
    
    % Plot's legends
    cb = colorbar;
    set(cb, 'Ticks', [0.25, 0.75], 'TickLabels', {'Desactivated (0)', 'Activated (1)'}, 'FontSize', 10);
    
    % function plotting
    subplot(1, 4, 4);
    plot(history_perf, 1:num_paliers_reals, '-ko', 'LineWidth', 1.5, 'MarkerFaceColor', [0 0 0], 'MarkerSize', 4);
    hold on;
    
    plot(best_perf_overall, num_paliers_reals, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [1 1 1]);
    plot(best_perf_overall, num_paliers_reals, 'rx', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    grid on; box on;
    set(gca, 'YDir', 'normal', 'TickDir', 'out', 'LineWidth', 1.2, 'YTickLabel', [], 'FontSize', 10);
    ylim([0.5, num_paliers_reals + 0.5]);
    xlabel('Performance P', 'FontSize', 12, 'FontWeight', 'bold');
    title('P(temp)', 'FontSize', 13, 'FontWeight', 'bold');
    
    hold off; 
    if other_figs == true
    figure('Name', 'Simulated Annealing - Convergence Diagnostics', 'Color', 'w', 'Position', [150, 150, 1000, 400]);
    tiledlayout(1, 3, 'TileSpacing', 'compact');
    
    % 1. Discrimanation performance
    nexttile
    plot(hist_iter_P, 'k', 'LineWidth', 1.2);
    hold on;
    plot(hist_iter_bestP, 'r', 'LineWidth', 2);
    xlabel('Iteration');
    ylabel('Performance P');
    legend({'Current', 'Best'}, 'Location', 'best');
    title('Discrimination performance');
    grid on; box on;
    
    % 2. Population size
    nexttile
    plot(hist_iter_size, 'b', 'LineWidth', 1.2);
    hold on;
    yline(length(best_subpop), '--r', 'Optimal found', 'LineWidth', 1.5);
    xlabel('Iteration');
    ylabel('Population size');
    title('Subpopulation Size Track');
    grid on; box on;
    
    % 3. Temperature cooling 
    nexttile
    semilogy(hist_iter_temp, 'm', 'LineWidth', 1.5);
    xlabel('Iteration');
    ylabel('Temperature');
    title('Cooling Schedule (Log Scale)');
    grid on; box on;
    
    % Global title
    sgtitle(sprintf('Simulated Annealing Dynamics | Best P = %.4f | Best Population = [%s]',...
        best_perf_overall, num2str(sort(best_subpop))));
    end
    shg;
end
end