%% Function to compute Labeled Line (LL) performance 
% Date: July 2026
% Author : Laure WOLFF 

function [Pn, Matrix] = calculate_labeled_line_P(CellMatrix, num_neurons, S, R, tmin, tmax, metric, plotting)
    
    num_trials = S * R;
    % Allocation of the 3D matrix (Trials x Trials x Neurons)
    Matrix = zeros(num_trials, num_trials, num_neurons);
    Pn = zeros(num_neurons, 1);
    
    for n = 1:num_neurons
                
        MatrixD = zeros(num_trials, num_trials);
        
        %% 1. Pre-extraction and ground-truth stimulus tracking
        Precomputed_Trains = cell(1, num_trials);
        vrai_stimulus = zeros(1, num_trials); % To record the true stimulus ID for each trial
        
        counter = 1;
        for st = 1:S
            for rp = 1:R
                if counter <= num_trials
                    % Extract the spike train
                    Precomputed_Trains{counter} = sort(CellMatrix{n, st, rp});
                    % Record the ground-truth stimulus ID for this trial
                    vrai_stimulus(counter) = st; 
                    counter = counter + 1;
                end
            end
        end
    
        %% 2. Pairwise Distance Calculation Matrix
        for t_a = 1:num_trials
            train_A = Precomputed_Trains{t_a};
            
            for t_b = (t_a + 1):num_trials 
                train_B = Precomputed_Trains{t_b};
                
                % SILENCE AND NOISE HANDLING WITHOUT ARTIFACTUAL RED LINES
                if isempty(train_A) && isempty(train_B)
                    dval = 0; % Two trials without spikes are identical
                elseif isempty(train_A) || isempty(train_B)
                    % If one train is empty, map to NaN for clean background masking
                    dval = NaN;
                else
                    % --- CLASSIC VECTORIZED COMPUTATION ---
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
                        dval = 0.5;
                    end
                end
    
                MatrixD(t_a, t_b) = dval;
                MatrixD(t_b, t_a) = dval;
            end
        end
        
        Matrix(:, :, n) = MatrixD;
    
        %% 3. Dynamic Global Performance Verification (P)
        profil_preference = zeros(1, num_trials);
        for t = 1:num_trials
            st_id = vrai_stimulus(t);
            if ~isempty(Precomputed_Trains{t})
                profil_preference(t) = st_id; 
            else
                profil_preference(t) = 0;     
            end
        end
        
        [Grid_A, Grid_B] = meshgrid(profil_preference, profil_preference);
        tri_upper_strict = triu(true(num_trials), 1);  
        
        is_intra = (Grid_A == Grid_B) & tri_upper_strict;
        is_inter = (Grid_A ~= Grid_B) & tri_upper_strict;
        
        if any(is_intra, 'all') && any(is_inter, 'all')
            % RECTIFICATION: Distance measures dissimilarity.
            % Intra-stimulus values are close to 0 (blue) and Inter-stimulus values are large (red/yellow).
            P = nanmean(MatrixD(is_inter)) - nanmean(MatrixD(is_intra)); 
        else
            P = 0;
        end
    
        if isnan(P), P = -Inf; end
        Pn(n) = P;
    end
    
    %% 4. PLOTTING SECTION
    if plotting == true
        trial_labels = cell(1, num_trials);
        counter = 1;
        for st = 1:S
            for rp = 1:R
                trial_labels{counter} = sprintf('S%d-R%d', st, rp);
                counter = counter + 1;
            end
        end
        
        figure('Name', 'LL Pairwise Distance Matrices Dn per Neuron', 'Color', 'w');
        
        cols = ceil(sqrt(num_neurons * 1.25)); 
        rows = ceil(num_neurons / cols);
        
        for n = 1:num_neurons
            subplot(rows, cols, n);
            
            h = imagesc(Matrix(:, :, n)); 
            colormap('jet'); 
            axis square;
            
            set(h, 'AlphaData', ~isnan(Matrix(:, :, n))); 
            set(gca, 'Color', [0.9 0.9 0.9]); 
            
            if num_neurons <= 6
                colorbar;
            end
            
            if num_trials <= 12
                set(gca, 'XTick', 1:num_trials, 'YTick', 1:num_trials);
                set(gca, 'XTickLabel', trial_labels, 'YTickLabel', trial_labels, 'FontSize', 6, 'TickLabelInterpreter', 'none');
                xtickangle(45);
            else
                set(gca, 'XTick', [1, num_trials], 'YTick', [1, num_trials]);
                set(gca, 'XTickLabel', {'1', num2str(num_trials)}, 'YTickLabel', {'1', num2str(num_trials)}, 'FontSize', 8);
            end
            
            title(sprintf('Neuron %d (P = %.2f)', n, Pn(n)), 'FontSize', 9, 'FontWeight', 'bold');
            xlabel('Trials', 'FontSize', 7); 
            ylabel('Trials', 'FontSize', 7);
        end
        shg;
    end
end