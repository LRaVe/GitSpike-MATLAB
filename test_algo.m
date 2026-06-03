%% Script Principal - Reproduction Figure 3 (Satuvuori 2018)
% Reprise à zéro - Approche Summed Population Strict
clear; clc; close all;

% Paramètres de la simulation
num_stimuli = 4;         % S
num_repetitions = 5;     % R
num_neurons = 7;         % N
num_coding_neurons = 3;  % Les 3 premiers codent, les autres sont du bruit
num_trials = num_stimuli * num_repetitions;

t1 = 0; t2 = 1;          
refrac = 0.01;           
base_rate_noise = 10;    % Bruit de fond (10 Hz)
base_rate_coding = 40;   % Signal (40 Hz)

metric_choice = 'ISI_ADAPTIVE'; 
rng(42); % Graine aléatoire fixe pour la reproductibilité

% 1. Génération des trains de spikes (Poisson)
CellMatrix = cell(num_neurons, num_stimuli, num_repetitions);

for st = 1:num_stimuli
    for rp = 1:num_repetitions
        for nc = 1:num_neurons
            
            if nc <= num_coding_neurons
                % Profil complémentaire : chaque neurone codant préfère certains stimuli
                if (nc == 1 && (st == 1 || st == 2)) || ...
                   (nc == 2 && (st == 2 || st == 3)) || ...
                   (nc == 3 && (st == 3 || st == 4))
                    rate = base_rate_coding;
                else
                    rate = base_rate_noise;
                end
            else
                % Neurones de bruit pur (10 Hz constants partout)
                rate = base_rate_noise;
            end
            
            approx_spikes = round((t2 - t1) * rate * 3) + 10;
            intervals = refrac - log(1 - rand(1, approx_spikes)) / rate;
            spikes = cumsum(intervals);
            CellMatrix{nc, st, rp} = spikes(spikes >= t1 & spikes <= t2);
        end
    end
end

% 2. Algorithme Bottom-Up (Summed Population)
subpop_choisie = zeros(1, num_neurons);     
subpop_dispo = 1:num_neurons;              
history_perf = zeros(1, num_neurons);

fprintf('Début de l''algorithme Bottom-Up (Summed Population)...\n');

for k = 1:num_neurons
    best_perf_etape = -Inf; 
    best_neurone_etape = -1;
    
    for i = 1:length(subpop_dispo)
        neurone_test = subpop_dispo(i);
        
        % Construction du masque de sélection [0 ou 1]
        selection = zeros(num_neurons, 1);
        if k > 1
            selection(subpop_choisie(1:k-1)) = 1;
        end
        selection(neurone_test) = 1;
        
        % Appel de la fonction de calcul (Fusion interne des spikes)
        [current_perf, ~] = calculate_integrated_P_optimized(...
            CellMatrix, selection, num_stimuli, num_repetitions, t1, t2, metric_choice);
        
        if current_perf > best_perf_etape
            best_perf_etape = current_perf;
            best_neurone_etape = neurone_test;
        end
    end
    
    subpop_choisie(k) = best_neurone_etape;
    subpop_dispo(subpop_dispo == best_neurone_etape) = [];
    history_perf(k) = best_perf_etape;
    
    fprintf('Taille %d : Neurone ajouté = %d | Performance Max P = %.4f\n', ...
        k, best_neurone_etape, best_perf_etape);
end

% 3. Affichage du graphique de droite (Figure 3)
figure('Name', 'Reproduction Figure 3', 'Position', [100, 100, 600, 450]);
plot(1:num_neurons, history_perf, '-or', 'LineWidth', 2, 'MarkerFaceColor', 'r');
grid on;
xlabel('Taille de la sous-population (k)');
ylabel('Performance de discrimination (P)');
title('Maximum de performance par taille (Bottom-Up SP)');