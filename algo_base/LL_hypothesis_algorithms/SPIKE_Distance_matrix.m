% %% Function to compute and plot pairwise SPIKE-distance matrices
% % Date: July 2026
% % Author: Laure WOLFF 
% 
% function All_Matrix_D = SPIKE_Distance_matrix(CellMatrix, num_neurons, S, R, tmin, tmax, metric, plotting)
% 
%     num_trials = S * R;
%     All_Matrix_D = zeros(num_trials, num_trials, num_neurons);
% 
% 
%     for n = 1:num_neurons
%         MatrixD = zeros(num_trials, num_trials);
% 
%         %% 1. Pre-extract spike trains and track stimulus identities
%         Precomputed_Trains = cell(1, num_trials);
%         vrai_stimulus = zeros(1, num_trials); 
% 
%         counter = 1;
%         for st = 1:S
%             for rp = 1:R
%                 if counter <= num_trials
%                     Precomputed_Trains{counter} = sort(CellMatrix{n, st, rp});
%                     vrai_stimulus(counter) = st; 
%                     counter = counter + 1;
%                 end
%             end
%         end
% 
%         % %% 2. Pairwise Distance Calculation Loop
%         % for t_a = 1:num_trials
%         %     train_A = Precomputed_Trains{t_a};
%         % 
%         %     for t_b = (t_a + 1):num_trials 
%         %         train_B = Precomputed_Trains{t_b};
%         % 
%         %         % Handle empty train edge cases (Silence vs Silence / Silence vs Signal)
%         %         if isempty(train_A) && isempty(train_B)
%         %             dval = 0; 
%         %         elseif isempty(train_A) || isempty(train_B)
%         %             %dval = NaN; % Marked as NaN for clean background masking
%         %             dval = 1; 
%         %         else
%         %             % --- SPIKE-DISTANCE Core Formulation ---
%         %             t_all = [tmin, sort(unique([train_A, train_B])), tmax];
%         %             t_diff = diff(t_all);
%         %             t_mids = (t_all(1:end-1) + t_all(2:end)) / 2;
%         % 
%         %             [~, ~, bin_A] = histcounts(t_mids, [-Inf, train_A, Inf]);
%         %             [~, ~, bin_B] = histcounts(t_mids, [-Inf, train_B, Inf]);
%         % 
%         %             if strcmpi(metric, 'SPIKE_DISTANCE')
%         %                 idx_p_A = max(1, bin_A - 1); idx_n_A = min(length(train_A), bin_A);
%         %                 x_p = train_A(idx_p_A); x_p(bin_A - 1 < 1) = tmin;
%         %                 x_a = train_A(idx_n_A); x_a(bin_A > length(train_A)) = tmax;
%         % 
%         %                 idx_p_B = max(1, bin_B - 1); idx_n_B = min(length(train_B), bin_B);
%         %                 y_p = train_B(idx_p_B); y_p(bin_B - 1 < 1) = tmin;
%         %                 y_a = train_B(idx_n_B); y_a(bin_B > length(train_B)) = tmax;
%         % 
%         %                 isi_x = x_a - x_p; 
%         %                 isi_y = y_a - y_p;
%         %                 dt_x_p = t_mids - x_p; 
%         %                 dt_x_a = x_a - t_mids;
%         %                 target_x = x_a; 
%         %                 target_x(dt_x_p < dt_x_a) = x_p(dt_x_p < dt_x_a);
%         % 
%         %                 dt_y_p = t_mids - y_p; 
%         %                 dt_y_a = y_a - t_mids;
%         %                 target_y = y_a; 
%         %                 target_y(dt_y_p < dt_y_a) = y_p(dt_y_p < dt_y_a);
%         % 
%         %                 [~, ~, bin_tX] = histcounts(target_x, [-Inf, train_B, Inf]);
%         %                 idx_p_Bt = max(1, bin_tX - 1); 
%         %                 idx_n_Bt = min(length(train_B), bin_tX);
%         %                 near_B = train_B(idx_n_Bt);
%         %                 use_p_B = (target_x - train_B(idx_p_Bt)) < (train_B(idx_n_Bt) - target_x) & (bin_tX - 1 >= 1);
%         %                 near_B(use_p_B) = train_B(idx_p_Bt(use_p_B));
%         %                 min_dxy = abs(target_x - near_B);
%         % 
%         %                 [~, ~, bin_tY] = histcounts(target_y, [-Inf, train_A, Inf]);
%         %                 idx_p_At = max(1, bin_tY - 1); 
%         %                 idx_n_At = min(length(train_A), bin_tY);
%         %                 near_A = train_A(idx_n_At);
%         %                 use_p_A = (target_y - train_A(idx_p_At)) < (train_A(idx_n_At) - target_y) & (bin_tY - 1 >= 1);
%         %                 near_A(use_p_A) = train_A(idx_p_At(use_p_A));
%         %                 min_dyx = abs(target_y - near_A);
%         % 
%         %                 S_x = (dt_x_p .* min_dyx + dt_x_a .* min_dyx) ./ isi_x;
%         %                 S_y = (dt_y_p .* min_dxy + dt_y_a .* min_dxy) ./ isi_y;
%         %                 S_t_list = (S_x .* isi_y + S_y .* isi_x) ./ ((isi_x + isi_y) .* max(isi_x, isi_y));
%         %                 dval = sum(S_t_list .* t_diff) / (tmax - tmin);
%         %             end
%         %         end
%         % 
%         %         MatrixD(t_a, t_b) = dval;
%         %         MatrixD(t_b, t_a) = dval;
%         %     end
%         % end
% 
%         %% 2. Pairwise Distance Calculation Loop
%         for t_a = 1:num_trials
%             train_A = Precomputed_Trains{t_a};
%             for t_b = (t_a + 1):num_trials 
%                 train_B = Precomputed_Trains{t_b};
% 
%                 % Handle empty train edge cases (Silence vs Silence / Silence vs Signal)
%                 if isempty(train_A) && isempty(train_B)
%                     dval = 0; 
%                 elseif isempty(train_A) || isempty(train_B)
%                     dval = 1; 
%                 else
%                     % --- SPIKE-DISTANCE Core Formulation ---
%                     if strcmpi(metric, 'SPIKE_DISTANCE')
% 
%                         % 1. Vecteur temps combiné (tous les points d'inflexion de la courbe)
%                         t_all = unique([tmin, train_A, train_B, tmax]);
%                         N_t = length(t_all);
% 
%                         % Initialisation des profils de distance temporelle
%                         S_x = zeros(1, N_t);
%                         S_y = zeros(1, N_t);
% 
%                         % 2. Calcul des profils pour chaque point de t_all
%                         for i = 1:N_t
%                             t = t_all(i);
% 
%                             % --- Train A ---
%                             % Trouver le spike précédent et suivant dans A
%                             idx_after_A = find(train_A >= t, 1, 'first');
%                             if isempty(idx_after_A)
%                                 x_p = train_A(end); x_a = tmax;
%                             elseif idx_after_A == 1
%                                 x_p = tmin; x_a = train_A(1);
%                             else
%                                 x_p = train_A(idx_after_A - 1); x_a = train_A(idx_after_A);
%                             end
%                             isi_x = x_a - x_p;
%                             dt_x_p = t - x_p;
%                             dt_x_a = x_a - t;
% 
%                             % Distance au spike le plus proche de A
%                             [min_val_A, ~] = min(abs(t - train_A));
% 
%                             % --- Train B ---
%                             % Trouver le spike précédent et suivant dans B
%                             idx_after_B = find(train_B >= t, 1, 'first');
%                             if isempty(idx_after_B)
%                                 y_p = train_B(end); y_a = tmax;
%                             elseif idx_after_B == 1
%                                 y_p = tmin; y_a = train_B(1);
%                             else
%                                 y_p = train_B(idx_after_B - 1); y_a = train_B(idx_after_B);
%                             end
%                             isi_y = y_a - y_p;
%                             dt_y_p = t - y_p;
%                             dt_y_a = y_a - t;
% 
%                             % Distance au spike le plus proche de B
%                             [min_val_B, ~] = min(abs(t - train_B));
% 
%                             % --- Calcul des distances locales de Kreuz S_x et S_y ---
%                             % S_x : distance de A par rapport à B
%                             S_x(i) = (dt_x_p * min_val_B + dt_x_a * min_val_B) / isi_x; % Note: se simplifie en min_val_B car (dt_x_p + dt_x_a) = isi_x
%                             % S_y : distance de B par rapport à A
%                             S_y(i) = (dt_y_p * min_val_A + dt_y_a * min_val_A) / isi_y; % Se simplifie en min_val_A
%                         end
% 
%                         % 3. Combinaison des profils S(t)
%                         % Calcul des ISI courants pour chaque point t_all pour la pondération
%                         % (On ré-obtient isi_x et isi_y pour tout le vecteur t_all)
%                         isi_x_all = zeros(1, N_t);
%                         isi_y_all = zeros(1, N_t);
%                         for i = 1:N_t
%                             t = t_all(i);
%                             % ISI x
%                             ia = find(train_A >= t, 1, 'first');
%                             if isempty(ia)
%                                 xp = train_A(end);
%                                 xa = tmax;
%                             elseif ia == 1
%                                 xp = tmin; 
%                                 xa = train_A(1);
%                             else
%                                 xp = train_A(ia-1);
%                                 xa = train_A(ia); 
%                             end
%                             isi_x_all(i) = xa - xp;
% 
%                             % ISI y
%                             ib = find(train_B >= t, 1, 'first');
%                             if isempty(ib)
%                                 yp = train_B(end); 
%                                 ya = tmax;
%                             elseif ib == 1
%                                 yp = tmin; 
%                                 ya = train_B(1);
%                             else
%                                 yp = train_B(ib-1); 
%                                 ya = train_B(ib); 
%                             end
%                             isi_y_all(i) = ya - yp;
%                         end
% 
%                         % Profil final S(t)
%                         S_t = (S_x .* isi_y_all + S_y .* isi_x_all) ./ ((isi_x_all + isi_y_all) .* max(isi_x_all, isi_y_all));
% 
%                         % 4. Intégration par la méthode des trapèzes
%                         dval = trapz(t_all, S_t) / (tmax - tmin);
%                     end
%                 end
%                 MatrixD(t_a, t_b) = dval;
%                 MatrixD(t_b, t_a) = dval;
%             end
%         end
% 
%         All_Matrix_D(:, :, n) = MatrixD;
%     end
% 
%     %% 3. PLOTTING SECTION
%     if plotting == true
%         trial_labels = cell(1, num_trials);
%         counter = 1;
%         for st = 1:S
%             for rp = 1:R
%                 trial_labels{counter} = sprintf('S%d-R%d', st, rp);
%                 counter = counter + 1;
%             end
%         end
% 
%         figure('Name', 'Pairwise SPIKE Distance Matrices per Neuron', 'Color', 'w');
% 
%         cols = ceil(sqrt(num_neurons * 1.25)); 
%         rows = ceil(num_neurons / cols);
% 
%         for n = 1:num_neurons
%             subplot(rows, cols, n);
% 
%             % Generate matrix visual profile
%             h = imagesc(All_Matrix_D(:, :, n)); 
%             colormap('jet'); 
%             axis square;
%             colorbar;
% 
%             % Background layout masking for NaNs (displays undefined silent pairings as light grey)
%             set(h, 'AlphaData', ~isnan(All_Matrix_D(:, :, n))); 
%             set(gca, 'Color', [0.9 0.9 0.9]); 
% 
%             if num_trials <= 12
%                 set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
%                 set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'FontSize', 6, 'TickLabelInterpreter', 'none');
%                 xtickangle(45);
%             else
%                 set(gca, 'XTick', [1, num_trials], 'YTick', [1, num_trials]);
%                 set(gca, 'XTickLabel', {'1', num2str(num_trials)}, 'YTickLabel', {'1', num2str(num_trials)}, 'FontSize', 8);
%             end
% 
%             title(sprintf('Neuron %d ', n), 'FontSize', 10, 'FontWeight', 'bold');
%             xlabel('Trials', 'FontSize', 8); 
%             ylabel('Trials', 'FontSize', 8);
%         end
%         shg;
%     end
% end
% 
% 
% 
% 
% 

%% Function to compute and plot pairwise SPIKE-distance matrices
% Date: July 2026
% Author: Laure WOLFF (Corrected Version)

function All_Matrix_D = SPIKE_Distance_matrix(CellMatrix, num_neurons, S, R, tmin, tmax, metric, plotting)

    num_trials = S * R;
    All_Matrix_D = zeros(num_trials, num_trials, num_neurons);

    for n = 1:num_neurons
        MatrixD = zeros(num_trials, num_trials);

        %% 1. Pre-extract spike trains
        Precomputed_Trains = cell(1, num_trials);
        counter = 1;
        for st = 1:S
            for rp = 1:R
                if counter <= num_trials
                    Precomputed_Trains{counter} = sort(CellMatrix{n, st, rp});
                    counter = counter + 1;
                end
            end
        end

        %% 2. Pairwise Distance Calculation Loop (Kreuz 2013 Exact Formulation)
        for t_a = 1:num_trials
            train_A = Precomputed_Trains{t_a};
            for t_b = (t_a + 1):num_trials 
                train_B = Precomputed_Trains{t_b};
                
                if isempty(train_A) && isempty(train_B)
                    dval = 0; 
                elseif isempty(train_A) || isempty(train_B)
                    dval = 1; 
                else
                    if strcmpi(metric, 'SPIKE_DISTANCE')
                        
                        % Grille temporelle combinée (points d'inflexion)
                        t_all = unique([tmin, train_A, train_B, tmax]);
                        N_t = length(t_all);
                        
                        S_x = zeros(1, N_t);
                        S_y = zeros(1, N_t);
                        isi_x_all = zeros(1, N_t);
                        isi_y_all = zeros(1, N_t);
                        
                        for i = 1:N_t
                            t = t_all(i);
                            
                            %% --- TRAIN A ---
                            idx_after_A = find(train_A >= t, 1, 'first');
                            if isempty(idx_after_A)
                                x_p = train_A(end); x_a = tmax;
                            %% Cas particulier des frontières (Bords)
                            elseif idx_after_A == 1
                                x_p = tmin; x_a = train_A(1);
                            else
                                x_p = train_A(idx_after_A - 1); x_a = train_A(idx_after_A);
                            end
                            isi_x = x_a - x_p;
                            isi_x_all(i) = isi_x;
                            
                            dt_x_p = t - x_p;
                            dt_x_a = x_a - t;
                            
                            %% --- TRAIN B ---
                            idx_after_B = find(train_B >= t, 1, 'first');
                            if isempty(idx_after_B)
                                y_p = train_B(end); y_a = tmax;
                            elseif idx_after_B == 1
                                y_p = tmin; y_a = train_B(1);
                            else
                                y_p = train_B(idx_after_B - 1); y_a = train_B(idx_after_B);
                            end
                            isi_y = y_a - y_p;
                            isi_y_all(i) = isi_y;
                            
                            dt_y_p = t - y_p;
                            dt_y_a = y_a - t;
                            
                            %% --- FORMULATION EXACTE DES DISTANCES AUX SPIKES VOISINS ---
                            % Distance du spike x_p au train B
                            d_xp_B = min(abs(x_p - train_B));
                            if isempty(d_xp_B), d_xp_B = 0; end
                            
                            % Distance du spike x_a au train B
                            d_xa_B = min(abs(x_a - train_B));
                            if isempty(d_xa_B), d_xa_B = 0; end
                            
                            % Distance du spike y_p au train A
                            d_yp_A = min(abs(y_p - train_A));
                            if isempty(d_yp_A), d_yp_A = 0; end
                            
                            % Distance du spike y_a au train A
                            d_ya_A = min(abs(y_a - train_A));
                            if isempty(d_ya_A), d_ya_A = 0; end
                            
                            %% --- PROFILS DIRECTIONNELS DE KREUZ ---
                            S_x(i) = (dt_x_p * d_xa_B + dt_x_a * d_xp_B) / isi_x; 
                            S_y(i) = (dt_y_p * d_ya_A + dt_y_a * d_yp_A) / isi_y; 
                        end
                        
                        % Combinaison bivariate finale pondérée par les ISIs
                        S_t = (S_x .* isi_y_all + S_y .* isi_x_all) ./ ((isi_x_all + isi_y_all) .* max(isi_x_all, isi_y_all));
                        
                        % Intégration trapézoïdale
                        dval = trapz(t_all, S_t) / (tmax - tmin);
                    end
                end
                MatrixD(t_a, t_b) = dval;
                MatrixD(t_b, t_a) = dval;
            end
        end
        All_Matrix_D(:, :, n) = MatrixD;
    end

    %% 3. PLOTTING SECTION
    if plotting == true
        trial_labels = cell(1, num_trials);
        counter = 1;
        for st = 1:S
            for rp = 1:R
                trial_labels{counter} = sprintf('S%d-R%d', st, rp);
                counter = counter + 1;
            end
        end

        figure('Name', 'Pairwise SPIKE Distance Matrices per Neuron', 'Color', 'w');
        cols = ceil(sqrt(num_neurons * 1.25)); 
        rows = ceil(num_neurons / cols);

        for n = 1:num_neurons
            subplot(rows, cols, n);
            h = imagesc(All_Matrix_D(:, :, n)); 
            colormap('jet'); 
            axis square;
            colorbar;
            
            % Ajustement dynamique de l'échelle de couleur [0, max]
            max_val = max(reshape(All_Matrix_D(:,:,n), [], 1));
            if max_val > 0
                clim([0 max_val]);
            else
                clim([0 1]);
            end

            set(h, 'AlphaData', ~isnan(All_Matrix_D(:, :, n))); 
            set(gca, 'Color', [0.9 0.9 0.9]); 

            if num_trials <= 12
                set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
                set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'FontSize', 6, 'TickLabelInterpreter', 'none');
                xtickangle(45);
            else
                set(gca, 'XTick', [1, num_trials], 'YTick', [1, num_trials]);
                set(gca, 'XTickLabel', {'1', num2str(num_trials)}, 'YTickLabel', {'1', num2str(num_trials)}, 'FontSize', 8);
            end

            title(sprintf('Neuron %d', n), 'FontSize', 10, 'FontWeight', 'bold');
            xlabel('Trials', 'FontSize', 8); 
            ylabel('Trials', 'FontSize', 8);
        end
        shg;
    end
end