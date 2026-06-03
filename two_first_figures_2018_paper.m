%% Main script to plot the two first figures of the 2018's paper
% Author: Laure WOLFF 
% Date: May 2026

clear; clc; close all;


% Initialization of the parameters
num_stimuli = 4;         % S : Number of stimuli
num_repetitions = 5;     % R : Number of the repetitions by stimuli 
num_neurons = 7;        % N : Global number of the neurons 
num_coding_neurons = 3;  % c : Number of the coding neurons 
num_trials = num_stimuli * num_repetitions;

% Dataset of the section 3 of the 2018 paper

t1 = 0; t2 = 1;          % Time windows
%refrac = 0.002;          % "absolute refractory period of 2 ms" 
refrac = 0.05;
base_rate_noise = 10;    
base_rate_coding = 30;  

metric_choice = 'ISI_ADAPTIVE'; 

% Data generation with the Poisson law (Pooled Coding Subpopulation)
CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);

    % Master Template Generation for each Stimulus (Pooled Response)
MasterTemplates = cell(num_coding_neurons, num_stimuli);
for st = 1:num_stimuli
    for nc = 1:num_coding_neurons
        rate = base_rate_coding + randi([-5, 5]); 
        approx_spikes = round((t2 - t1) * rate * 3) + 10;
        intervals = local_f_poisson(approx_spikes, rate, refrac);
        spikes = cumsum(intervals);
        MasterTemplates{nc, st} = spikes(spikes >= t1 & spikes <= t2);
    end
end

    % Filling the global trial cell matrix
for st = 1:num_stimuli
    for rp = 1:num_repetitions
        for nc = 1:num_neurons
            if nc <= num_coding_neurons
                CellMatrix{nc, st, rp} = MasterTemplates{nc, st};
            else
                approx_spikes = round((t2 - t1) * base_rate_noise * 3) + 10;
                intervals = local_f_poisson(approx_spikes, base_rate_noise, refrac);
                spikes = cumsum(intervals);
                CellMatrix{nc, st, rp} = spikes(spikes >= t1 & spikes <= t2);
            end
        end
    end
end

% Creation of the labels
trial_labels = cell(1, num_trials);
counter = 1;
for st = 1:num_stimuli
    for rp = 1:num_repetitions
        trial_labels{counter} = sprintf('S%d-R%d', st, rp);
        counter = counter + 1;
    end
end

% Plotting the global raster plot
figure('Name', 'Global Dataset Raster Plot', 'Position', [100, 100, 950, 650]);
hold on;

color_noise = [0.0 0.45 0.74];  % Blue for the non-coding neurons (NC)
color_coding  = [0.85 0.33 0.1];  % Red for the coding neurons (C)

for t_idx = 1:num_trials
    st = floor((t_idx-1)/num_repetitions) + 1;
    rp = mod((t_idx-1), num_repetitions) + 1;
    
    for nc = 1:num_neurons
        spikes = CellMatrix{nc, st, rp};
        if ~isempty(spikes)
            if nc <= num_coding_neurons
                current_color = color_coding; line_width = 1.5;
            else
                current_color = color_noise; line_width = 1.0;
            end
            for sp = 1:length(spikes)
                line([spikes(sp), spikes(sp)], [t_idx - 0.4, t_idx + 0.4], ...
                     'Color', current_color, 'LineWidth', line_width);
            end
        end
    end
    line([t1, t2], [t_idx, t_idx], 'Color', [0.94 0.94 0.94], 'LineWidth', 0.5, 'LineStyle', '-');
end

box on; grid on;
set(gca, 'XGrid', 'on', 'YGrid', 'off');
xlim([t1, t2]); ylim([0.5, num_trials + 0.5]);
set(gca, 'YTick', 1:num_trials, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none', 'FontSize', 9);
xlabel('Time', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Trials (Stimuli & Repetitions)', 'FontSize', 11, 'FontWeight', 'bold');
title('Global Spike Train Raster Plot (Red: Coding | Blue: Noise)', 'FontSize', 13, 'FontWeight', 'bold');

for st = 1:(num_stimuli-1)
    sep_line = st * num_repetitions + 0.5;
    line([t1, t2], [sep_line, sep_line], 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2, 'LineStyle', '--');
end

% Creation of the three matrix 
coding_selection = [ones(num_coding_neurons, 1); zeros(num_neurons - num_coding_neurons, 1)];
noise_selection  = [zeros(num_coding_neurons, 1); ones(num_neurons - num_coding_neurons, 1)];
full_selection   = ones(num_neurons, 1);

%fprintf('--- Calcul optimisé des 3 Matrices (%d x %d essais) ---\n', num_trials, num_trials);
%tic;
[perf_A, Matrix_A] = calculate_integrated_P_optimized(CellMatrix, coding_selection, num_stimuli, num_repetitions, t1, t2, metric_choice);
[perf_B, Matrix_B] = calculate_integrated_P_optimized(CellMatrix, noise_selection, num_stimuli, num_repetitions, t1, t2, metric_choice);
[perf_C, Matrix_C] = calculate_integrated_P_optimized(CellMatrix, full_selection, num_stimuli, num_repetitions, t1, t2, metric_choice);
%toc;

% Plotting the three distance matrix 
figure('Name', 'SP Distances matrix', 'Position', [25, 150, 1500, 450]);

% Matrix A : Coding
subplot(1, 3, 1);
imagesc(Matrix_A); colormap('jet'); colorbar; axis square;
set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none');
xtickangle(45);
title(sprintf('A. Coding subpopulation (C)\nP = %.4f', perf_A), 'Color', 'r', 'FontWeight', 'bold');
xlabel('Trials'); ylabel('Trials');

% Matrix B : Non-Coding 
subplot(1, 3, 2);
imagesc(Matrix_B); colormap('jet'); colorbar; axis square;
set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none');
xtickangle(45);
title(sprintf('B. Non-coding subpopulation (NC)\nP = %.4f', perf_B), 'Color', 'b', 'FontWeight', 'bold');
xlabel('Trials');

% Matrix C : All
subplot(1, 3, 3);
imagesc(Matrix_C); colormap('jet'); colorbar; axis square;
set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none');
xtickangle(45);
title(sprintf('C. Full population (All)\nP = %.4f', perf_C), 'Color', 'k', 'FontWeight', 'bold');
xlabel('Trials');

hold off; shg;


%% =========================================================================
%% LOCAL FUNCTIONS BLOCK
%% =========================================================================
function poiss = local_f_poisson(len, rate, refrac)
    uniform = rand(1, len);
    poiss = refrac - log(1 - uniform) / rate;
end

function [P, MatrixD] = calculate_integrated_P_optimized(CellMatrix, selection, S, R, tmin, tmax, metric)
    idx_selected = find(selection == 1);
    num_trials = S * R;
    MatrixD = zeros(num_trials, num_trials);
    
    if isempty(idx_selected), P = -Inf; return; end
    
    % Optimized loop exploiting the symmetry of the matrix
    for t_a = 1:num_trials
        st_a = floor((t_a-1)/R) + 1;
        rp_a = mod(t_a-1, R) + 1;
        
        for t_b = (t_a + 1):num_trials 
            st_b = floor((t_b-1)/R) + 1;
            rp_b = mod(t_b-1, R) + 1;
            
            % Summed population hypothesis
            train_A = sort([CellMatrix{idx_selected, st_a, rp_a}]);
            train_B = sort([CellMatrix{idx_selected, st_b, rp_b}]);
            
            % A-ISI distance (or the ISI if MRTS = 0)
            if strcmpi(metric, 'ISI_ADAPTIVE')
                t_all = [tmin, sort([train_A, train_B]), tmax]; 
                It_list = zeros(1, length(t_all)-1);
                
                sum_sqr = 0; n_isi = 0;
                if length(train_A) >= 2
                    sum_sqr = sum_sqr + sum(diff(train_A).^2);
                    n_isi = n_isi + length(train_A) - 1;
                end
                if length(train_B) >= 2
                    sum_sqr = sum_sqr + sum(diff(train_B).^2);
                    n_isi = n_isi + length(train_B) - 1;
                end
                MRTS = 0; 
                if n_isi > 0, MRTS = (sum_sqr/n_isi)^0.5; end
                
                for k = 1 : length(t_all)-1
                    t_mid = (t_all(k) + t_all(k+1)) / 2;
                    
                    % Train A
                    if isempty(train_A)
                        vx = tmax - tmin;
                    elseif t_mid < train_A(1)
                        vx = train_A(1) - tmin;
                    elseif t_mid > train_A(end)
                        vx = tmax - train_A(end);
                    else
                        idx_v = find(train_A <= t_mid, 1, 'last');
                        if idx_v == length(train_A)
                            vx = tmax - train_A(end);
                        else
                            vx = train_A(idx_v+1) - train_A(idx_v); 
                        end
                    end
                    
                    % Train B
                    if isempty(train_B)
                        vy = tmax - tmin;
                    elseif t_mid < train_B(1)
                        vy = train_B(1) - tmin;
                    elseif t_mid > train_B(end)
                        vy = tmax - train_B(end);
                    else
                        idy = find(train_B <= t_mid, 1, 'last'); 
                        if idy == length(train_B) 
                            vy = tmax - train_B(end);
                        else
                            vy = train_B(idy+1) - train_B(idy); 
                        end
                    end
                    
                    It_list(k) = abs(vx - vy) / max([vx, vy, MRTS]);
                end
                dval = sum(It_list .* diff(t_all)) / (tmax - tmin);
            else
                dval = 0.5;
            end
            
            % Using the symetry of the matrix to assign
            MatrixD(t_a, t_b) = dval;
            MatrixD(t_b, t_a) = dval;

        end
    end
    
    % Final evaluation of global performance P
    sum_intra = 0; count_intra = 0;
    sum_inter = 0; count_inter = 0;
    for t_a = 1:num_trials
        st_a = floor((t_a-1)/R) + 1;
        for t_b = (t_a+1):num_trials
            st_b = floor((t_b-1)/R) + 1;
            dist_val = MatrixD(t_a, t_b);
            if st_a == st_b
                sum_intra = sum_intra + dist_val;
                count_intra = count_intra + 1;
            else
                sum_inter = sum_inter + dist_val;
                count_inter = count_inter + 1;
            end
        end
    end
    P = (sum_inter / count_inter) - (sum_intra / count_intra);
end