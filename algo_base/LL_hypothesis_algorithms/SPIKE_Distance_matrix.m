%% Function to compute and plot pairwise SPIKE-distance matrices
% Date: July 2026
% Author: Laure WOLFF 

function All_Matrix_D = SPIKE_Distance_matrix(CellMatrix, num_neurons, S, R, tmin, tmax, metric, showing, plotting)
    
    num_trials = S * R;
    % Preallocate the 3D distance matrix (Trials x Trials x Neurons)
    All_Matrix_D = zeros(num_trials, num_trials, num_neurons);
    
   
    for n = 1:num_neurons
        MatrixD = zeros(num_trials, num_trials);
        
        %% 1. Pre-extract spike trains and track stimulus identities
        Precomputed_Trains = cell(1, num_trials);
        vrai_stimulus = zeros(1, num_trials); 
        
        counter = 1;
        for st = 1:S
            for rp = 1:R
                if counter <= num_trials
                    Precomputed_Trains{counter} = sort(CellMatrix{n, st, rp});
                    vrai_stimulus(counter) = st; 
                    counter = counter + 1;
                end
            end
        end
    
        %% 2. Pairwise Distance Calculation Loop
        for t_a = 1:num_trials
            train_A = Precomputed_Trains{t_a};
            
            for t_b = (t_a + 1):num_trials 
                train_B = Precomputed_Trains{t_b};
                
                % Handle empty train edge cases (Silence vs Silence / Silence vs Signal)
                if isempty(train_A) && isempty(train_B)
                    dval = 0; 
                elseif isempty(train_A) || isempty(train_B)
                    dval = NaN; % Marked as NaN for clean background masking
                else
                    % --- SPIKE-DISTANCE Core Formulation ---
                    t_all = [tmin, sort(unique([train_A, train_B])), tmax];
                    t_diff = diff(t_all);
                    t_mids = (t_all(1:end-1) + t_all(2:end)) / 2;
    
                    [~, ~, bin_A] = histcounts(t_mids, [-Inf, train_A, Inf]);
                    [~, ~, bin_B] = histcounts(t_mids, [-Inf, train_B, Inf]);
    
                    if strcmpi(metric, 'SPIKE_DISTANCE')
                        idx_p_A = max(1, bin_A - 1); idx_n_A = min(length(train_A), bin_A);
                        x_p = train_A(idx_p_A); x_p(bin_A - 1 < 1) = tmin;
                        x_a = train_A(idx_n_A); x_a(bin_A > length(train_A)) = tmax;
                        
                        idx_p_B = max(1, bin_B - 1); idx_n_B = min(length(train_B), bin_B);
                        y_p = train_B(idx_p_B); y_p(bin_B - 1 < 1) = tmin;
                        y_a = train_B(idx_n_B); y_a(bin_B > length(train_B)) = tmax;
                        
                        isi_x = x_a - x_p; isi_y = y_a - y_p;
                        dt_x_p = t_mids - x_p; dt_x_a = x_a - t_mids;
                        target_x = x_a; target_x(dt_x_p < dt_x_a) = x_p(dt_x_p < dt_x_a);
                        
                        dt_y_p = t_mids - y_p; dt_y_a = y_a - t_mids;
                        target_y = y_a; target_y(dt_y_p < dt_y_a) = y_p(dt_y_p < dt_y_a);
                        
                        [~, ~, bin_tX] = histcounts(target_x, [-Inf, train_B, Inf]);
                        idx_p_Bt = max(1, bin_tX - 1); idx_n_Bt = min(length(train_B), bin_tX);
                        near_B = train_B(idx_n_Bt);
                        use_p_B = (target_x - train_B(idx_p_Bt)) < (train_B(idx_n_Bt) - target_x) & (bin_tX - 1 >= 1);
                        near_B(use_p_B) = train_B(idx_p_Bt(use_p_B));
                        min_dxy = abs(target_x - near_B);
                        
                        [~, ~, bin_tY] = histcounts(target_y, [-Inf, train_A, Inf]);
                        idx_p_At = max(1, bin_tY - 1); idx_n_At = min(length(train_A), bin_tY);
                        near_A = train_A(idx_n_At);
                        use_p_A = (target_y - train_A(idx_p_At)) < (train_A(idx_n_At) - target_y) & (bin_tY - 1 >= 1);
                        near_A(use_p_A) = train_A(idx_p_At(use_p_A));
                        min_dyx = abs(target_y - near_A);
                        
                        S_x = (dt_x_p .* min_dyx + dt_x_a .* min_dyx) ./ isi_x;
                        S_y = (dt_y_p .* min_dxy + dt_y_a .* min_dxy) ./ isi_y;
                        S_t_list = (S_x .* isi_y + S_y .* isi_x) ./ ((isi_x + isi_y) .* max(isi_x, isi_y));
                        dval = sum(S_t_list .* t_diff) / (tmax - tmin);
                        
                    elseif strcmpi(metric, 'ISI_ADAPTIVE')
                        sum_sqr = 0; n_isi = 0;
                        if length(train_A) >= 2, sum_sqr = sum_sqr + sum(diff(train_A).^2); n_isi = n_isi + length(train_A) - 1; end
                        if length(train_B) >= 2, sum_sqr = sum_sqr + sum(diff(train_B).^2); n_isi = n_isi + length(train_B) - 1; end
                        MRTS = 0; if n_isi > 0, MRTS = sqrt(sum_sqr/n_isi); end
                        
                        idx_v = bin_A - 1; vx = zeros(1, length(t_mids)); vx(idx_v < 1) = train_A(1) - tmin; vx(idx_v >= length(train_A)) = tmax - train_A(end);
                        valid_x = (idx_v >= 1) & (idx_v < length(train_A)); if any(valid_x), vx(valid_x) = train_A(idx_v(valid_x) + 1) - train_A(idx_v(valid_x)); end
                        
                        idy = bin_B - 1; vy = zeros(1, length(t_mids)); vy(idy < 1) = train_B(1) - tmin; vy(idy >= length(train_B)) = tmax - train_B(end);
                        valid_y = (idy >= 1) & (idy < length(train_B)); if any(valid_y), vy(valid_y) = train_B(idy(valid_y) + 1) - train_B(idy(valid_y)); end
                        
                        max_v = max([vx; vy; repmat(MRTS, 1, length(t_mids))], [], 1);
                        It_list = abs(vx - vy) / max_v;
                        dval = sum(It_list .* t_diff) / (tmax - tmin);
                    else
                        dval = 0.5; % Default fallback value
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
        
        figure('Name', 'Pairwise SPIKE Distance Matrices per Neuron', 'Color', 'w', 'Position', [100, 100, 1200, 800]);
        
        cols = ceil(sqrt(num_neurons * 1.25)); 
        rows = ceil(num_neurons / cols);
        
        for n = 1:num_neurons
            subplot(rows, cols, n);
            
            % Generate matrix visual profile
            h = imagesc(All_Matrix_D(:, :, n)); 
            colormap('jet'); 
            axis square;
            colorbar;
            
            % Background layout masking for NaNs (displays undefined silent pairings as light grey)
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
            
            title(sprintf('Neuron %d ', n), 'FontSize', 10, 'FontWeight', 'bold');
            xlabel('Trials', 'FontSize', 8); 
            ylabel('Trials', 'FontSize', 8);
        end
        shg;
    end
end