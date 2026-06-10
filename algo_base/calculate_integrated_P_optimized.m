%% Script to compute the performance in summed population hypothesis 
% Date: May-June 2026
% Author : Laure WOLFF

function [P, MatrixD] = calculate_integrated_P_optimized(CellMatrix, selection, S, R, tmin, tmax, metric)

    idx_selected = find(selection == 1);
    num_trials = S * R;
    MatrixD = zeros(num_trials, num_trials);
   
    if isempty(idx_selected), P = -Inf; return; end

    %% 1. Pré-extraction des trains pour la sous-population sélectionnée
    % On fusionne (Summed Population) pour chaque essai à l'avance pour éviter de le refaire dans la boucle imbriquée
    Precomputed_Trains = cell(1, num_trials);

    for t = 1:num_trials
        st = floor((t-1)/R) + 1;
        rp = mod(t-1, R) + 1;
        Precomputed_Trains{t} = sort([CellMatrix{idx_selected, st, rp}]);
    end

    %% 2. Distance's calcul
    for t_a = 1:num_trials
        train_A = Precomputed_Trains{t_a};
        is_empty_A = isempty(train_A);
       
        for t_b = (t_a + 1):num_trials 
            train_B = Precomputed_Trains{t_b};
            is_empty_B = isempty(train_B);
           
            if is_empty_A && is_empty_B
                dval = 0;
            elseif is_empty_A || is_empty_B
                dval = 1;
            else

                %% Spike distance (trying to avoid a loops to optimized the complexity) 
                if strcmpi(metric, 'SPIKE_DISTANCE')
                    t_all = [tmin, sort(unique([train_A, train_B])), tmax];
                    t_diff = diff(t_all);
                    t_mids = (t_all(1:end-1) + t_all(2:end)) / 2;
                    
                    % 1. Find previous and next spikes for ALL t_mids simultaneously
                    [~, ~, bin_A] = histcounts(t_mids, [-Inf, train_A, Inf]);
                    [~, ~, bin_B] = histcounts(t_mids, [-Inf, train_B, Inf]);
                    
                    % Bound and extract indices for train A (handling boundary conditions)
                    idx_p_A = bin_A - 1; idx_p_A(idx_p_A < 1) = 1;
                    idx_n_A = bin_A;     idx_n_A(idx_n_A > length(train_A)) = length(train_A);
                    x_p = train_A(idx_p_A); x_p(bin_A - 1 < 1) = tmin;
                    x_a = train_A(idx_n_A); x_a(bin_A > length(train_A)) = tmax;
                    
                    % Bound and extract indices for train B (handling boundary conditions)
                    idx_p_B = bin_B - 1; idx_p_B(idx_p_B < 1) = 1;
                    idx_n_B = bin_B;     idx_n_B(idx_n_B > length(train_B)) = length(train_B);
                    y_p = train_B(idx_p_B); y_p(bin_B - 1 < 1) = tmin;
                    y_a = train_B(idx_n_B); y_a(bin_B > length(train_B)) = tmax;
                    
                    % Calculate the ISI distance
                    isi_x = x_a - x_p;
                    isi_y = y_a - y_p;
                    
                    % 2. Determine target_x and target_y for the entire vector
                    dt_x_p = t_mids - x_p;
                    dt_x_a = x_a - t_mids;
                    target_x = x_a; 
                    target_x(dt_x_p < dt_x_a) = x_p(dt_x_p < dt_x_a);
                    
                    dt_y_p = t_mids - y_p;
                    dt_y_a = y_a - t_mids;
                    target_y = y_a; 
                    target_y(dt_y_p < dt_y_a) = y_p(dt_y_p < dt_y_a);
                    
                    % 3. Nearest neighbor interpolation
                    % Find the closest spike in train_B to target_x using vectorized extrapolation
                    % (On utilise une astuce d'indexation pour simuler le plus proche voisin de manière vectorisée)
                    min_dxy = abs(target_x - interp1(train_B, train_B, target_x, 'nearest', 'extrap'));
                    min_dyx = abs(target_y - interp1(train_A, train_A, target_y, 'nearest', 'extrap'));
                    
                    % 4. Compute vectorized distance profiles
                    S_x = (dt_x_p .* min_dyx + dt_x_a .* min_dyx) ./ isi_x;
                    S_y = (dt_y_p .* min_dxy + dt_y_a .* min_dxy) ./ isi_y;
                    S_t_list = (S_x .* isi_y + S_y .* isi_x) ./ ((isi_x + isi_y) .* max(isi_x, isi_y));
                    
                    % Time integration 
                    dval = sum(S_t_list .* t_diff) / (tmax - tmin);


                %% ISI distance
                elseif strcmpi(metric, 'ISI_ADAPTIVE')
                    t_all = [tmin, sort([train_A, train_B]), tmax]; 
                    t_diff = diff(t_all);
                    len_t = length(t_all) - 1;
                    It_list = zeros(1, len_t);
                    t_mids = (t_all(1:end-1) + t_all(2:end)) / 2;              
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

                    if n_isi > 0 
                        MRTS = sqrt(sum_sqr/n_isi); 
                    end
                    
                    [~, ~, bin_A] = histcounts(t_mids, [-Inf, train_A, Inf]);
                    [~, ~, bin_B] = histcounts(t_mids, [-Inf, train_B, Inf]);
                
                    for k = 1 : len_t
                        % Train A
                        idx_v = bin_A(k) - 1;
                        if idx_v < 1
                            vx = train_A(1) - tmin;
                        elseif idx_v >= length(train_A)
                            vx = tmax - train_A(end);
                        else
                            vx = train_A(idx_v+1) - train_A(idx_v); 
                        end
                       
                        % Train B
                        idy = bin_B(k) - 1;
                        if idy < 1
                            vy = train_B(1) - tmin;
                        elseif idy >= length(train_B)
                            vy = tmax - train_B(end);
                        else
                            vy = train_B(idy+1) - train_B(idy); 
                        end
                       
                        It_list(k) = abs(vx - vy) / max([vx, vy, MRTS]);
                    end

                    dval = sum(It_list .* t_diff) / (tmax - tmin);
                else
                    dval = 0.5;
                end
            end
           
            MatrixD(t_a, t_b) = dval;
            MatrixD(t_b, t_a) = dval;
        end
    end
 
    %% 3. Calcul of the performance P
    stim_A = floor(((1:num_trials)-1)/R) + 1;
    [Grid_A, Grid_B] = meshgrid(stim_A, stim_A);
    
    tri_upper = triu(true(num_trials), 1);  
    is_intra = (Grid_A == Grid_B) & tri_upper;
    is_inter = (Grid_A ~= Grid_B) & tri_upper;

    P = mean(MatrixD(is_inter)) - mean(MatrixD(is_intra));
end