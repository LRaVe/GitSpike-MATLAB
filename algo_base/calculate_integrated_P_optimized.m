%% Script to compute the performance in summed population hypothesis
% Date: May-June 2026
% Author : Laure WOLFF

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
            
            %% Summed Population hypothesis
            train_A = sort([CellMatrix{idx_selected, st_a, rp_a}]);
            train_B = sort([CellMatrix{idx_selected, st_b, rp_b}]);
            
            if strcmpi(metric, 'SPIKE_DISTANCE')
                t_all = [tmin, sort(unique([train_A, train_B])), tmax];
                S_t_list = zeros(1, length(t_all)-1);
                
                if isempty(train_A) && isempty(train_B)
                    dval = 0;
                elseif isempty(train_A) || isempty(train_B)
                    dval = 1;
                else
                    for k = 1 : length(t_all)-1
                        t_mid = (t_all(k) + t_all(k+1)) / 2;
                        
                        % Train A 
                        idx_a = find(train_A > t_mid, 1, 'first');
                        idx_p = find(train_A <= t_mid, 1, 'last');
                        
                        if isempty(idx_p), x_p = tmin; else, x_p = train_A(idx_p); end
                        if isempty(idx_a), x_a = tmax; else, x_a = train_A(idx_a); end
                        isi_x = x_a - x_p;
                        
                        % time distance to the A spike
                        dt_x_p = t_mid - x_p;
                        dt_x_a = x_a - t_mid;
                        
                        % Train B 
                        idy_a = find(train_B > t_mid, 1, 'first');
                        idy_p = find(train_B <= t_mid, 1, 'last');
                        
                        if isempty(idy_p), y_p = tmin; else, y_p = train_B(idy_p); end
                        if isempty(idy_a), y_a = tmax; else, y_a = train_B(idy_a); end
                        isi_y = y_a - y_p;
                        
                        % time distance to the B spike
                        dt_y_p = t_mid - y_p;
                        dt_y_a = y_a - t_mid;
                        
                        % For the nearest peak in A, find its distance to the nearest peak in B
                        [~, closest_A_idx] = min([dt_x_p, dt_x_a]);
                        if closest_A_idx == 1, target_x = x_p; else, target_x = x_a; end
                        [~, id_b_closest] = min(abs(train_B - target_x));
                        min_dxy = abs(target_x - train_B(id_b_closest));
                        
                        % For the nearest peak in B, find its distance to
                        % the nearest peak in A
                        [~, closest_B_idx] = min([dt_y_p, dt_y_a]);
                        if closest_B_idx == 1, target_y = y_p; else, target_y = y_a; end
                        [~, id_a_closest] = min(abs(train_A - target_y));
                        min_dyx = abs(target_y - train_A(id_a_closest));
                        
                        % Calcul time profile S_x and S_y
                        S_x = (dt_x_p * min_dyx + dt_x_a * min_dyx) / isi_x; 
                        S_y = (dt_y_p * min_dxy + dt_y_a * min_dxy) / isi_y;
                        
                        % Local inter-spike interval (ISI) weighting
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