function [best_subpop, max_global_perf, history] = f_bottom_up_discrimination(N, num_stimuli, distance_matrix_cell, num_repetitions)
    % Algorithme Bottom-Up basé sur la moyenne des matrices de distance pré-calculées
    
    subpop_choisie = zeros(1, N);     
    subpop_dispo = 1:N;              
    
    history_perf = zeros(1, N);
    history_subpops = cell(1, N);
    
    fprintf('Début de l''algorithme Bottom-Up...\n');
    
    for k = 1:N
        best_perf_etape = -Inf; 
        best_neurone_etape = -1;
        
        for i = 1:length(subpop_dispo)
            neurone_test = subpop_dispo(i);
            
            % Liste temporaire des neurones testés à cette étape
            subpop_temp = [subpop_choisie(1:k-1), neurone_test];
            
            % Combinaison par moyenne des matrices individuelles
            Matrix_combo = zeros(size(distance_matrix_cell{1}));
            for n = 1:length(subpop_temp)
                Matrix_combo = Matrix_combo + distance_matrix_cell{subpop_temp(n)};
            end
            Matrix_combo = Matrix_combo / length(subpop_temp);
            
            % Évaluation locale et ultra-rapide du score P
            current_perf = local_evaluer_P(Matrix_combo, num_stimuli, num_repetitions);
            
            if current_perf > best_perf_etape
                best_perf_etape = current_perf;
                best_neurone_etape = neurone_test;
            end
        end
        
        % Validation du meilleur neurone pour l'étape k
        subpop_choisie(k) = best_neurone_etape;
        subpop_dispo(subpop_dispo == best_neurone_etape) = [];
        
        history_perf(k) = best_perf_etape;
        history_subpops{k} = subpop_choisie(1:k);
        
        fprintf('Taille %d : Neurone ajouté = %d | Performance Max P = %.4f\n', k, best_neurone_etape, best_perf_etape);
    end
    
    [max_global_perf, idx_best_taille] = max(history_perf);
    best_subpop = history_subpops{idx_best_taille};
    
    history.perf = history_perf;
    history.subpops = history_subpops;
    
    fprintf('\n=== RÉSULTAT FINAL BOTTOM-UP ===\n');
    fprintf('La sous-population idéale a une taille de %d neurones.\n', idx_best_taille);
    fprintf('Indices des neurones codants détectés : %s\n', num2str(best_subpop));
    fprintf('Performance maximale P : %.4f\n', max_global_perf);
end

%% =========================================================================
%% FONCTION LOCALE : Calcul de P à partir d'une matrice de distance combinée
%% =========================================================================
function P = local_evaluer_P(MatrixD, S, R)
    num_trials = S * R;
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