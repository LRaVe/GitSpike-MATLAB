function [P, MatrixD] = calculate_integrated_P_optimized(CellMatrix, selection, S, R, tmin, tmax, metric)
    idx_selected = find(selection == 1);
    num_trials = S * R;
    MatrixD = zeros(num_trials, num_trials);
    
    if isempty(idx_selected), P = -Inf; return; end
    
    % Boucle optimisée exploitant la symétrie de la matrice
    for t_a = 1:num_trials
        st_a = floor((t_a-1)/R) + 1;
        rp_a = mod(t_a-1, R) + 1;
        
        for t_b = (t_a + 1):num_trials 
            st_b = floor((t_b-1)/R) + 1;
            rp_b = mod(t_b-1, R) + 1;
            
            total_dist = 0;
            for n = 1:length(idx_selected)
                neuron_idx = idx_selected(n);
                train_A = CellMatrix{neuron_idx, st_a, rp_a};
                train_B = CellMatrix{neuron_idx, st_b, rp_b};
                
                % Calcul de la distance ISI Adaptative
                if strcmpi(metric, 'ISI_ADAPTIVE')
                    t_all = [tmin, sort([train_A, train_B]), tmax]; 
                    It_list = zeros(1, length(t_all)-1);
                    
                    sum_sqr = 0; 
                    n_isi = 0;
                    
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
                        MRTS = (sum_sqr/n_isi)^0.5; 
                    end
                    
                    for k = 1 : length(t_all)-1
                        t_mid = (t_all(k) + t_all(k+1)) / 2;
                        
                        % Train A
                        if isempty(train_A)
                            vx = tmax - tmin;
                        else
                            idx_v = find(train_A <= t_mid, 1, 'last');
                            if isempty(idx_v)
                                vx = train_A(1) - tmin;
                            elseif idx_v == length(train_A)
                                vx = tmax - train_A(end);
                            else
                                vx = train_A(idx_v+1) - train_A(idx_v); 
                            end
                        end
                        
                        % Train B
                        if isempty(train_B)
                            vy = tmax - tmin;
                        else
                            idy = find(train_B <= t_mid, 1, 'last'); 
                            if isempty(idy)
                                vy = train_B(1) - tmin;
                            elseif idy == length(train_B) 
                                vy = tmax - train_B(end);
                            else
                                vy = train_B(idy+1) - train_B(idy); 
                            end
                        end
                        
                        It_list(k) = abs(vx - vy) / max([vx, vy, MRTS]);
                    end
                    dval = sum(It_list .* diff(t_all)) / (tmax - tmin);
                else
                    dval = 0.5;
                end
                total_dist = total_dist + dval;
            end
            
            % Attribution symétrique (Summed Population)
            val_final = total_dist / length(idx_selected);
            MatrixD(t_a, t_b) = val_final;
            MatrixD(t_b, t_a) = val_final; 
        end
    end
    
    % Évaluation finale de la performance globale P
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