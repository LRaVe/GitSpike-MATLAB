%% Main script to plot the two first figures of the 2018's paper
% Author: Laure WOLFF 
% Date: May 2026

clear; clc; close all;

% Initialisation des paramètres
num_stimuli = 4;         % S : Nombre de stimuli
num_repetitions = 5;     % R : Nombre de répétitions par stimulus 
num_neurons = 7;         % N : Nombre global de neurones 
num_coding_neurons = 3;  % c : Nombre de neurones codants 
num_trials = num_stimuli * num_repetitions;

t1 = 0; t2 = 1;          % Fenêtre temporelle
refrac = 0.05;           % Période réfractaire
base_rate_noise = 10;    
base_rate_coding = 30;  
jitter_std = 0.01;       % Écart-type du bruit temporel entre répétitions (10 ms)

% --- Génération du Dataset (Pooled Coding Subpopulation) ---
CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);
MasterTemplates = cell(num_coding_neurons, num_stimuli);

% 1. Génération des Master Templates pour chaque Stimulus
for st = 1:num_stimuli
    for nc = 1:num_coding_neurons
        rate = base_rate_coding + randi([-5, 5]); 
        approx_spikes = round((t2 - t1) * rate * 3) + 10;
        intervals = local_f_poisson(approx_spikes, rate, refrac);
        spikes = cumsum(intervals);
        MasterTemplates{nc, st} = spikes(spikes >= t1 & spikes <= t2);
    end
end

% 2. Remplissage de la matrice globale des essais (avec Jitter stochastique)
for st = 1:num_stimuli
    for rp = 1:num_repetitions
        for nc = 1:num_neurons
            if nc <= num_coding_neurons
                % REPRODUCTION CONFORME : On ajoute une gigue temporelle (bruit) 
                % propre à chaque répétition pour éviter des trains 100% identiques.
                template_spikes = MasterTemplates{nc, st};
                jittered_spikes = template_spikes + jitter_std * randn(size(template_spikes));
                % On s'assure que les pics restent dans les bornes t1 et t2
                jittered_spikes = sort(jittered_spikes(jittered_spikes >= t1 & jittered_spikes <= t2));
                CellMatrix{nc, st, rp} = jittered_spikes;
            else
                % Neurones non-codants (Bruit pur indépendant à chaque essai)
                approx_spikes = round((t2 - t1) * base_rate_noise * 3) + 10;
                intervals = local_f_poisson(approx_spikes, base_rate_noise, refrac);
                spikes = cumsum(intervals);
                CellMatrix{nc, st, rp} = spikes(spikes >= t1 & spikes <= t2);
            end
        end
    end
end

% Création des étiquettes d'essais
trial_labels = cell(1, num_trials);
counter = 1;
for st = 1:num_stimuli
    for rp = 1:num_repetitions
        trial_labels{counter} = sprintf('S%d-R%d', st, rp);
        counter = counter + 1;
    end
end

%% =========================================================================
%% BLOC D'AFFICHAGE DES RASTER PLOTS INDÉPENDANTS (FENÊTRE UNIQUE)
%% =========================================================================
color_coding  = [0.85 0.33 0.1];  % Orange/Rouge pour les codants
color_noise   = [0.0 0.45 0.74];  % Bleu pour le bruit

% Configuration commune des sous-graphiques
plot_types = {'CODING', 'NON-CODING', 'FULL'};
titles = {'1. Coding Subpopulation Only (Neurons 1 to 3)', ...
          '2. Non-Coding Subpopulation Only (Neurons 4 to 7)', ...
          '3. Full Population (All Neurons)'};

% Création d'une seule fenêtre verticale [gauche, bas, largeur, hauteur]
figure('Name', 'Subpopulations Spike Train Raster Plots', 'Position', [200, 50, 900, 900]);

for f = 1:3
    % Crée une grille de 3 lignes et 1 colonne, et active la ligne 'f'
    subplot(3, 1, f);
    hold on;
    
    for t_idx = 1:num_trials
        st = floor((t_idx-1)/num_repetitions) + 1;
        rp = mod((t_idx-1), num_repetitions) + 1;
        
        % Détermination des neurones à tracer selon le sous-graphique
        if strcmp(plot_types{f}, 'CODING')
            neurons_to_plot = 1:num_coding_neurons;
        elseif strcmp(plot_types{f}, 'NON-CODING')
            neurons_to_plot = (num_coding_neurons + 1):num_neurons;
        else
            neurons_to_plot = 1:num_neurons;
        end
        
        % Tracé des pics pour les neurones sélectionnés
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
        % Ligne de base pour chaque essai
        line([t1, t2], [t_idx, t_idx], 'Color', [0.94 0.94 0.94], 'LineWidth', 0.5);
    end
    
    % Mise en forme esthétique du graphique
    box on; grid on;
    set(gca, 'XGrid', 'on', 'YGrid', 'off');
    xlim([t1, t2]); ylim([0.5, num_trials + 0.5]);
    
    % Optionnel : On n'affiche les labels de l'axe Y (essais) de manière détaillée 
    % que si nécessaire, pour ne pas surcharger visuellement la figure.
    set(gca, 'YTick', 1:num_trials, 'YTickLabel', trial_labels, 'TickLabelInterpreter', 'none', 'FontSize', 8);
    
    % On n'affiche le label X complet que sur le graphique du bas pour épurer le rendu
    if f == 3
        xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
    end
    
    ylabel('Trials', 'FontSize', 10, 'FontWeight', 'bold');
    title(titles{f}, 'FontSize', 11, 'FontWeight', 'bold');
    
    % Lignes de démarcation entre les différents stimuli (S1, S2, S3, S4)
    for st_sep = 1:(num_stimuli-1)
        sep_line = st_sep * num_repetitions + 0.5;
        line([t1, t2], [sep_line, sep_line], 'Color', [0.1 0.1 0.1], 'LineWidth', 1.2, 'LineStyle', '--');
    end
    hold off;
end

shg;

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
    
    for t_a = 1:num_trials
        st_a = floor((t_a-1)/R) + 1;
        rp_a = mod(t_a-1, R) + 1;
        
        for t_b = (t_a + 1):num_trials 
            st_b = floor((t_b-1)/R) + 1;
            rp_b = mod(t_b-1, R) + 1;
            
            % Summed Population hypothesis
            train_A = sort([CellMatrix{idx_selected, st_a, rp_a}]);
            train_B = sort([CellMatrix{idx_selected, st_b, rp_b}]);
            
            if strcmpi(metric, 'SPIKE_DISTANCE')
                t_all = [tmin, sort(unique([train_A, train_B])), tmax];
                S_t_list = zeros(1, length(t_all)-1);
                
                % Si l'un des deux trains est complètement vide
                if isempty(train_A) && isempty(train_B)
                    dval = 0;
                elseif isempty(train_A) || isempty(train_B)
                    dval = 1; % Distance maximale si un train est vide
                else
                    for k = 1 : length(t_all)-1
                        t_mid = (t_all(k) + t_all(k+1)) / 2;
                        
                        % TRAIN A 
                        idx_a = find(train_A > t_mid, 1, 'first');
                        idx_p = find(train_A <= t_mid, 1, 'last');
                        
                        if isempty(idx_p), x_p = tmin; else, x_p = train_A(idx_p); end
                        if isempty(idx_a), x_a = tmax; else, x_a = train_A(idx_a); end
                        isi_x = x_a - x_p;
                        
                        % Distances temporelles instantanées aux pics de A
                        dt_x_p = t_mid - x_p;
                        dt_x_a = x_a - t_mid;
                        
                        % TRAIN B 
                        idy_a = find(train_B > t_mid, 1, 'first');
                        idy_p = find(train_B <= t_mid, 1, 'last');
                        
                        if isempty(idy_p), y_p = tmin; else, y_p = train_B(idy_p); end
                        if isempty(idy_a), y_a = tmax; else, y_a = train_B(idy_a); end
                        isi_y = y_a - y_p;
                        
                        % Distances temporelles instantanées aux pics de B
                        dt_y_p = t_mid - y_p;
                        dt_y_a = y_a - t_mid;
                        
                        % --- CO-LOCALISATION ( Kreuz et al. ) ---
                        % Pour le pic le plus proche de A, trouver sa distance au pic le plus proche dans B
                        [~, closest_A_idx] = min([dt_x_p, dt_x_a]);
                        if closest_A_idx == 1, target_x = x_p; else, target_x = x_a; end
                        [~, id_b_closest] = min(abs(train_B - target_x));
                        min_dxy = abs(target_x - train_B(id_b_closest));
                        
                        % Pour le pic le plus proche de B, trouver sa distance au pic le plus proche dans A
                        [~, closest_B_idx] = min([dt_y_p, dt_y_a]);
                        if closest_B_idx == 1, target_y = y_p; else, target_y = y_a; end
                        [~, id_a_closest] = min(abs(train_A - target_y));
                        min_dyx = abs(target_y - train_A(id_a_closest));
                        
                        % --- CALCUL DES PROFILS TEMPORELS S_x ET S_y ---
                        S_x = (dt_x_p * min_dyx + dt_x_a * min_dyx) / isi_x; 
                        S_y = (dt_y_p * min_dxy + dt_y_a * min_dxy) / isi_y;
                        
                        % --- COMBINAISON ET NORMALISATION ---
                        % Pondération par l'inter-pic local (ISI)
                        S_t_list(k) = (S_x * isi_y + S_y * isi_x) / ((isi_x + isi_y) * max([isi_x, isi_y]));
                    end
                    % Time integration
                    dval = sum(S_t_list .* diff(t_all)) / (tmax - tmin);
                end

            elseif strcmpi(metric, 'ISI_ADAPTIVE')
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
                        if idx_v == length(train_A), vx = tmax - train_A(end);
                        else, vx = train_A(idx_v+1) - train_A(idx_v); end
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
                        if idy == length(train_B), vy = tmax - train_B(end);
                        else, vy = train_B(idy+1) - train_B(idy); end
                    end
                    It_list(k) = abs(vx - vy) / max([vx, vy, MRTS]);
                end
                dval = sum(It_list .* diff(t_all)) / (tmax - tmin);
            else
                dval = 0.5;
            end
            
            % Using the symetrie to fill the matrix
            MatrixD(t_a, t_b) = dval;
            MatrixD(t_b, t_a) = dval;
        end
    end
    
    % Calculate the global performace P
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