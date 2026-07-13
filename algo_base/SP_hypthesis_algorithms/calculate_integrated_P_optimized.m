%% Script to compute the performance in summed population hypothesis 
% Date: May-June 2026
% Author : Laure WOLFF

function [P, MatrixD] = calculate_integrated_P_optimized(CellMatrix, selection, S, R, tmin, tmax, metric)
    idx_selected = find(selection == 1);
    num_trials = S * R;
    MatrixD = zeros(num_trials, num_trials);

    if isempty(idx_selected), P = -Inf; return; end

    %% 1. Pre-extraction and summation of spike trains per trial
    Precomputed_Trains = cell(1, num_trials);
    for t = 1:num_trials
        st = floor((t-1)/R) + 1;
        rp = mod(t-1, R) + 1;
        % Population Pooling (Summed Population)
        Precomputed_Trains{t} = sort([CellMatrix{idx_selected, st, rp}]);
    end

    %% 2. Pairwise Distance Calculation Matrix
    for t_a = 1:num_trials
        train_A = Precomputed_Trains{t_a};
        is_empty_A = isempty(train_A);
        for t_b = (t_a + 1):num_trials 
            train_B = Precomputed_Trains{t_b};
            is_empty_B = isempty(train_B);
            
            % Handle empty train edge cases (Silence vs Silence / Silence vs Signal)
            if is_empty_A && is_empty_B
                dval = 0;
            elseif is_empty_A || is_empty_B
                dval = 1;
            else
                % Combined time vector containing all unique inflexion points (spikes + boundaries)
                t_all = [tmin, sort(unique([train_A, train_B])), tmax];
                N_t = length(t_all);
                t_diff = diff(t_all); % Width of each constant segment
                
                %% =============================================================
                %% OPTION 1: SPIKE DISTANCE (Exact Continuous Integration)
                %% =============================================================
                if strcmpi(metric, 'SPIKE_DISTANCE')
                    % Initialize instantaneous distance profiles
                    S_x = zeros(1, N_t);
                    S_y = zeros(1, N_t);
                    
                    % Evaluate distances at each critical time point to handle piecewise linearity
                    for i = 1:N_t
                        t = t_all(i);
                        
                        % --- Train A Profiling ---
                        % Find enclosing spikes in Train A
                        idx_after_A = find(train_A >= t, 1, 'first');
                        if isempty(idx_after_A)
                            x_p = train_A(end); x_a = tmax;
                        elseif idx_after_A == 1
                            x_p = tmin; x_a = train_A(1);
                        else
                            x_p = train_A(idx_after_A - 1); x_a = train_A(idx_after_A);
                        end
                        isi_x = x_a - x_p;
                        dt_x_p = t - x_p;
                        dt_x_a = x_a - t;
                        
                        % Distance to the nearest spike in Train A
                        min_val_A = min(abs(t - train_A));
                        
                        % --- Train B Profiling ---
                        % Find enclosing spikes in Train B
                        idx_after_B = find(train_B >= t, 1, 'first');
                        if isempty(idx_after_B)
                            y_p = train_B(end); y_a = tmax;
                        elseif idx_after_B == 1
                            y_p = tmin; y_a = train_B(1);
                        else
                            y_p = train_B(idx_after_B - 1); y_a = train_B(idx_after_B);
                        end
                        isi_y = y_a - y_p;
                        dt_y_p = t - y_p;
                        dt_y_a = y_a - t;
                        
                        % Distance to the nearest spike in Train B
                        min_val_B = min(abs(t - train_B));
                        
                        % --- Local Kreuz Distances ---
                        S_x(i) = (dt_x_p * min_val_B + dt_x_a * min_val_B) / isi_x; 
                        S_y(i) = (dt_y_p * min_val_A + dt_y_a * min_val_A) / isi_y; 
                    end
                    
                    % Recompute instantaneous ISI vectors over t_all for bivariate weighting
                    isi_x_all = zeros(1, N_t);
                    isi_y_all = zeros(1, N_t);
                    for i = 1:N_t
                        t = t_all(i);
                        ia = find(train_A >= t, 1, 'first');
                        if isempty(ia)
                            xp = train_A(end); 
                            xa = tmax;
                        elseif ia == 1
                            xp = tmin; 
                            xa = train_A(1);
                        else
                            xp = train_A(ia-1); 
                            xa = train_A(ia); 
                        end
                        isi_x_all(i) = xa - xp;
                        
                        ib = find(train_B >= t, 1, 'first');
                        if isempty(ib)
                            yp = train_B(end); 
                            ya = tmax;
                        elseif ib == 1
                            yp = tmin; 
                            ya = train_B(1);
                        else
                            yp = train_B(ib-1); 
                            ya = train_B(ib); 
                        end
                        isi_y_all(i) = ya - yp;
                    end
                    
                    % Calculate final bivariate profile S(t)
                    S_t = (S_x .* isi_y_all + S_y .* isi_x_all) ./ ((isi_x_all + isi_y_all) .* max(isi_x_all, isi_y_all));
                    
                    % Integration using the trapezoidal method
                    dval = trapz(t_all, S_t) / (tmax - tmin);
                    
                %% =============================================================
                %% OPTION 2: ISI ADAPTIVE DISTANCE (Vectorized Step-wise Integration)
                %% =============================================================
                elseif strcmpi(metric, 'ISI_ADAPTIVE')
                    % Compute Mean Absolute Interval Deviation (MRTS)
                    sum_sqr = 0; n_isi = 0;
                    if length(train_A) >= 2
                        sum_sqr = sum_sqr + sum(diff(train_A).^2);
                        n_isi = n_isi + length(train_A) - 1;
                    end
                    if length(train_B) >= 2
                        sum_sqr = sum_sqr + sum(diff(train_B).^2);
                        n_isi = n_isi + length(train_B) - 1;
                    end
                    MRTS = 0; if n_isi > 0, MRTS = sqrt(sum_sqr/n_isi); end
                    
                    % Extract ISI values based on the start of each geometric segment.
                    % Using t_all(1:end-1) guarantees correct staircase step assignment.
                    [~, ~, bin_A] = histcounts(t_all(1:end-1), [-Inf, train_A, Inf]);
                    [~, ~, bin_B] = histcounts(t_all(1:end-1), [-Inf, train_B, Inf]);
                    
                    % Vectorized instantaneous ISI extraction for Train A
                    idx_v = bin_A - 1;
                    vx = zeros(1, length(t_diff));
                    vx(idx_v < 1) = train_A(1) - tmin;
                    vx(idx_v >= length(train_A)) = tmax - train_A(end);
                    valid_x = (idx_v >= 1) & (idx_v < length(train_A));
                    if any(valid_x)
                        vx(valid_x) = train_A(idx_v(valid_x) + 1) - train_A(idx_v(valid_x));
                    end
                    
                    % Vectorized instantaneous ISI extraction for Train B
                    idy = bin_B - 1;
                    vy = zeros(1, length(t_diff));
                    vy(idy < 1) = train_B(1) - tmin;
                    vy(idy >= length(train_B)) = tmax - train_B(end);
                    valid_y = (idy >= 1) & (idy < length(train_B));
                    if any(valid_y)
                        vy(valid_y) = train_B(idy(valid_y) + 1) - train_B(idy(valid_y));
                    end
                    
                    % Compute the adaptive profile and execute direct rectangular summation.
                    % Since ISIs are perfectly constant across each t_diff segment, this is mathematically exact.
                    max_v = max([vx; vy; repmat(MRTS, 1, length(t_diff))], [], 1);
                    It_list = abs(vx - vy) / max_v;
                    
                    dval = sum(It_list .* t_diff) / (tmax - tmin);
                else
                    dval = 0.5;
                end
            end
            % Fill symmetric distance matrix
            MatrixD(t_a, t_b) = dval;
            MatrixD(t_b, t_a) = dval;
        end
    end

    %% 3. Mathematical Global Performance Verification (P)
    stim_A = floor(((1:num_trials)-1)/R) + 1;
    [Grid_A, Grid_B] = meshgrid(stim_A, stim_A);

    tri_upper = triu(true(num_trials), 1);  
    is_intra = (Grid_A == Grid_B) & tri_upper;
    is_inter = (Grid_A ~= Grid_B) & tri_upper;

    P = mean(MatrixD(is_inter)) - mean(MatrixD(is_intra));

    % Safeguard against extreme mathematical noise or empty outputs
    if isnan(P), P = -Inf; end
end


% %% Script to compute the performance in summed population hypothesis 
% % Date: June 2026
% % Author : Laure WOLFF
% function [P, MatrixD] = calculate_integrated_P_optimized(CellMatrix, selection, S, R, tmin, tmax, metric)
%     idx_selected = find(selection == 1);
%     num_trials = S * R;
%     MatrixD = zeros(num_trials, num_trials);
% 
%     if isempty(idx_selected), P = -Inf; return; end
% 
%     %% 1. Pre-extraction and summation of spike trains per trial
%     Precomputed_Trains = cell(1, num_trials);
%     for t = 1:num_trials
%         st = floor((t-1)/R) + 1;
%         rp = mod(t-1, R) + 1;
%         Precomputed_Trains{t} = sort([CellMatrix{idx_selected, st, rp}]);
%     end
% 
%     %% 2. Pairwise Distance Calculation Matrix
%     for t_a = 1:num_trials
%         train_A = Precomputed_Trains{t_a};
%         is_empty_A = isempty(train_A);
% 
%         for t_b = (t_a + 1):num_trials 
%             train_B = Precomputed_Trains{t_b};
%             is_empty_B = isempty(train_B);
% 
%             if is_empty_A && is_empty_B
%                 dval = 0;
%             elseif is_empty_A || is_empty_B
%                 dval = 1;
%             else
%                 %% =============================================================
%                 %% OPTION 1: SPIKE DISTANCE (Votre cœur vectorisé autonome)
%                 %% =============================================================
%                 if strcmpi(metric, 'SPIKE_DISTANCE')
%                     % Grille temporelle combinée et nettoyée du bruit numérique
%                     t_all = unique([tmin, train_A, train_B, tmax]);
%                     % Sécurité : éliminer les points virtuellement identiques pour éviter isi = 0
%                     t_all = t_all([true, diff(t_all) > 1e-9]);
% 
%                     if t_all(end) < tmax
%                         t_all = [t_all, tmax]; 
%                     end
% 
%                     % Localisation des intervalles avec histcounts standard
%                     [~, ~, bin_A] = histcounts(t_all, [-Inf, train_A, Inf]);
%                     [~, ~, bin_B] = histcounts(t_all, [-Inf, train_B, Inf]);
% 
%                     % Intervalles et bornes pour le train A
%                     idx_p_A = max(1, bin_A - 1); idx_n_A = min(length(train_A), bin_A);
%                     x_p = train_A(idx_p_A); x_p(bin_A - 1 < 1) = tmin;
%                     x_a = train_A(idx_n_A); x_a(bin_A > length(train_A)) = tmax;
%                     isi_x = x_a - x_p;
%                     isi_x(isi_x <= 0) = 1e-9; % Protection stricte anti-NaN
% 
%                     % Intervalles et bornes pour le train B
%                     idx_p_B = max(1, bin_B - 1); idx_n_B = min(length(train_B), bin_B);
%                     y_p = train_B(idx_p_B); y_p(bin_B - 1 < 1) = tmin;
%                     y_a = train_B(idx_n_B); y_a(bin_B > length(train_B)) = tmax;
%                     isi_y = y_a - y_p;
%                     isi_y(isi_y <= 0) = 1e-9; % Protection stricte anti-NaN
% 
%                     % Distances temporelles instantanées aux spikes encadrants
%                     dt_x_p = t_all - x_p; dt_x_a = x_a - t_all;
%                     dt_y_p = t_all - y_p; dt_y_a = y_a - t_all;
% 
%                     % Recherche vectorisée du plus proche voisin depuis chaque point de t_all
%                     try
%                         idx_near_B = knnsearch(train_B', t_all');
%                         min_dxy = abs(t_all - train_B(idx_near_B));
%                     catch
%                         min_dxy = zeros(size(t_all));
%                         for i = 1:length(t_all)
%                             [~, min_idx] = min(abs(t_all(i) - train_B));
%                             min_dxy(i) = abs(t_all(i) - train_B(min_idx));
%                         end
%                     end
% 
%                     try
%                         idx_near_A = knnsearch(train_A', t_all');
%                         min_dyx = abs(t_all - train_A(idx_near_A));
%                     catch
%                         min_dyx = zeros(size(t_all));
%                         for i = 1:length(t_all)
%                             [~, min_idx] = min(abs(t_all(i) - train_A));
%                             min_dyx(i) = abs(t_all(i) - train_A(min_idx));
%                         end
%                     end
% 
%                     % Profils de distances symétriques
%                     S_1 = (min_dyx .* dt_x_a + min_dyx .* dt_x_p) ./ isi_x;
%                     S_2 = (min_dxy .* dt_y_a + min_dxy .* dt_y_p) ./ isi_y;
% 
%                     % Formule exacte de la SPIKE-distance classique
%                     denom = (isi_x + isi_y).^2;
%                     S_t_list = (2 * (S_1 .* isi_y + S_2 .* isi_x)) ./ denom;
% 
%                     % Intégration par la méthode des trapèzes sur la vraie géométrie du signal
%                     dval = trapz(t_all, S_t_list) / (tmax - tmin);
% 
%                 %% =============================================================
%                 %% OPTION 2: ISI ADAPTIVE DISTANCE
%                 %% =============================================================
%                 elseif strcmpi(metric, 'ISI_ADAPTIVE')
%                     t_all = [tmin, sort(unique([train_A, train_B])), tmax];
%                     t_diff = diff(t_all);
%                     t_mids = (t_all(1:end-1) + t_all(2:end)) / 2;
% 
%                     [~, ~, bin_A] = histcounts(t_mids, [-Inf, train_A, Inf]);
%                     [~, ~, bin_B] = histcounts(t_mids, [-Inf, train_B, Inf]);
% 
%                     sum_sqr = 0; n_isi = 0;
%                     if length(train_A) >= 2, sum_sqr = sum_sqr + sum(diff(train_A).^2); n_isi = n_isi + length(train_A) - 1; end
%                     if length(train_B) >= 2, sum_sqr = sum_sqr + sum(diff(train_B).^2); n_isi = n_isi + length(train_B) - 1; end
%                     MRTS = 0; if n_isi > 0, MRTS = sqrt(sum_sqr/n_isi); end
% 
%                     idx_v = bin_A - 1;
%                     vx = zeros(1, length(t_mids));
%                     vx(idx_v < 1) = train_A(1) - tmin;
%                     vx(idx_v >= length(train_A)) = tmax - train_A(end);
%                     valid_x = (idx_v >= 1) & (idx_v < length(train_A));
%                     if any(valid_x), vx(valid_x) = train_A(idx_v(valid_x) + 1) - train_A(idx_v(valid_x)); end
% 
%                     idy = bin_B - 1;
%                     vy = zeros(1, length(t_mids));
%                     vy(idy < 1) = train_B(1) - tmin;
%                     vy(idy >= length(train_B)) = tmax - train_B(end);
%                     valid_y = (idy >= 1) & (idy < length(train_B));
%                     if any(valid_y), vy(valid_y) = train_B(idy(valid_y) + 1) - train_B(idy(valid_y)); end
% 
%                     max_v = max([vx; vy; repmat(MRTS, 1, length(t_mids))], [], 1);
%                     It_list = abs(vx - vy) / max_v;
% 
%                     dval = sum(It_list .* t_diff) / (tmax - tmin);
%                 else
%                     dval = 0.5;
%                 end
%             end
% 
%             MatrixD(t_a, t_b) = dval;
%             MatrixD(t_b, t_a) = dval;
%         end
%     end
% 
%     %% 3. Global Performance Verification (P)
%     stim_A = floor(((1:num_trials)-1)/R) + 1;
%     [Grid_A, Grid_B] = meshgrid(stim_A, stim_A);
% 
%     tri_upper = triu(true(num_trials), 1);  
%     is_intra = (Grid_A == Grid_B) & tri_upper;
%     is_inter = (Grid_A ~= Grid_B) & tri_upper;
% 
%     P = mean(MatrixD(is_inter)) - mean(MatrixD(is_intra));
%     if isnan(P), P = -Inf; end
% end