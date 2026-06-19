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
           
            if is_empty_A && is_empty_B
                dval = 0;
            elseif is_empty_A || is_empty_B
                dval = 1;
            else
                t_all = [tmin, sort(unique([train_A, train_B])), tmax];
                t_diff = diff(t_all);
                t_mids = (t_all(1:end-1) + t_all(2:end)) / 2;
                
                % Shared pre-computation: find enclosing spikes for all interval midpoints
                [~, ~, bin_A] = histcounts(t_mids, [-Inf, train_A, Inf]);
                [~, ~, bin_B] = histcounts(t_mids, [-Inf, train_B, Inf]);
                
                %% =============================================================
                %% OPTION 1: SPIKE DISTANCE (Fully Vectorized via histcounts)
                %% =============================================================
                if strcmpi(metric, 'SPIKE_DISTANCE')
                    % Handling boundaries for Train A
                    idx_p_A = max(1, bin_A - 1);
                    idx_n_A = min(length(train_A), bin_A);
                    x_p = train_A(idx_p_A); x_p(bin_A - 1 < 1) = tmin;
                    x_a = train_A(idx_n_A); x_a(bin_A > length(train_A)) = tmax;
                    
                    % Handling boundaries for Train B
                    idx_p_B = max(1, bin_B - 1);
                    idx_n_B = min(length(train_B), bin_B);
                    y_p = train_B(idx_p_B); y_p(bin_B - 1 < 1) = tmin;
                    y_a = train_B(idx_n_B); y_a(bin_B > length(train_B)) = tmax;
                    
                    isi_x = x_a - x_p;
                    isi_y = y_a - y_p;
                    
                    % Locate target spikes (closest spike in same train)
                    dt_x_p = t_mids - x_p; dt_x_a = x_a - t_mids;
                    target_x = x_a; target_x(dt_x_p < dt_x_a) = x_p(dt_x_p < dt_x_a);
                    
                    dt_y_p = t_mids - y_p; dt_y_a = y_a - t_mids;
                    target_y = y_a; target_y(dt_y_p < dt_y_a) = y_p(dt_y_p < dt_y_a);
                    
                    % Vectorized Nearest Neighbor Search for min_dxy (Target X inside Train B)
                    [~, ~, bin_tX] = histcounts(target_x, [-Inf, train_B, Inf]);
                    idx_p_Bt = max(1, bin_tX - 1); idx_n_Bt = min(length(train_B), bin_tX);
                    near_B = train_B(idx_n_Bt);
                    use_p_B = (target_x - train_B(idx_p_Bt)) < (train_B(idx_n_Bt) - target_x) & (bin_tX - 1 >= 1);
                    near_B(use_p_B) = train_B(idx_p_Bt(use_p_B));
                    min_dxy = abs(target_x - near_B);
                    
                    % Vectorized Nearest Neighbor Search for min_dyx (Target Y inside Train A)
                    [~, ~, bin_tY] = histcounts(target_y, [-Inf, train_A, Inf]);
                    idx_p_At = max(1, bin_tY - 1); idx_n_At = min(length(train_A), bin_tY);
                    near_A = train_A(idx_n_At);
                    use_p_A = (target_y - train_A(idx_p_At)) < (train_A(idx_n_At) - target_y) & (bin_tY - 1 >= 1);
                    near_A(use_p_A) = train_A(idx_p_At(use_p_A));
                    min_dyx = abs(target_y - near_A);
                    
                    % Final profiles and integration
                    S_x = (dt_x_p .* min_dyx + dt_x_a .* min_dyx) ./ isi_x;
                    S_y = (dt_y_p .* min_dxy + dt_y_a .* min_dxy) ./ isi_y;
                    S_t_list = (S_x .* isi_y + S_y .* isi_x) ./ ((isi_x + isi_y) .* max(isi_x, isi_y));
                    
                    dval = sum(S_t_list .* t_diff) / (tmax - tmin);
                    
                %% =============================================================
                %% OPTION 2: ISI ADAPTIVE DISTANCE (Fully Vectorized)
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
                    
                    % Vectorized instantaneous ISI extraction for Train A
                    idx_v = bin_A - 1;
                    vx = zeros(1, length(t_mids));
                    vx(idx_v < 1) = train_A(1) - tmin;
                    vx(idx_v >= length(train_A)) = tmax - train_A(end);
                    valid_x = (idx_v >= 1) & (idx_v < length(train_A));
                    if any(valid_x)
                        vx(valid_x) = train_A(idx_v(valid_x) + 1) - train_A(idx_v(valid_x));
                    end
                    
                    % Vectorized instantaneous ISI extraction for Train B
                    idy = bin_B - 1;
                    vy = zeros(1, length(t_mids));
                    vy(idy < 1) = train_B(1) - tmin;
                    vy(idy >= length(train_B)) = tmax - train_B(end);
                    valid_y = (idy >= 1) & (idy < length(train_B));
                    if any(valid_y)
                        vy(valid_y) = train_B(idy(valid_y) + 1) - train_B(idy(valid_y));
                    end
                    
                    % Calculation of the adaptive profile without loops
                    max_v = max([vx; vy; repmat(MRTS, 1, length(t_mids))], [], 1);
                    It_list = abs(vx - vy) / max_v;
                    
                    dval = sum(It_list .* t_diff) / (tmax - tmin);
                else
                    dval = 0.5;
                end
            end
           
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