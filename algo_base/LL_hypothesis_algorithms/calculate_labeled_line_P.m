%% Script to compute the performance in Labeled Line (LL) hypothesis
% Date: June 2026
% Author : Laure WOLFF
function [P_LL, MatrixD_LL] = calculate_labeled_line_P(CellMatrix, selection, S, R, tmin, tmax)
% CALCULATE_LABELED_LINE_P Compute performance under the Labeled Line hypothesis
% This function averages the distance matrices computed individually for each neuron.
    idx_selected = find(selection == 1);
    num_trials = S * R;
    num_neurons = length(idx_selected);
    
    % Preallocate the final output distance matrix
    MatrixD_LL = zeros(num_trials, num_trials);
   
    if isempty(idx_selected)
        P_LL = -Inf; 
        return; 
    end
    
    %% 1. Pairwise Distance Calculation Matrix (True LL Algorithm)
    for t_a = 1:num_trials
        st_a = floor((t_a-1)/R) + 1; rp_a = mod(t_a-1, R) + 1;
       
        for t_b = (t_a + 1):num_trials 
            st_b = floor((t_b-1)/R) + 1; rp_b = mod(t_b-1, R) + 1;
            
            % --- STEP 1: Create a master time grid pooling all spikes for this trial pair ---
            % 1.1 Compute the total number of spikes across all selected neurons for both trials
            total_spikes = 2; % Start with 2 slots for tmin and tmax
            for n = 1:num_neurons
                total_spikes = total_spikes + length(CellMatrix{idx_selected(n), st_a, rp_a}) ...
                                            + length(CellMatrix{idx_selected(n), st_b, rp_b});
            end
            
            % 1.2 Preallocate the master vector with fixed size
            all_spikes_pair = zeros(1, total_spikes);
            all_spikes_pair(1) = tmin;
            all_spikes_pair(2) = tmax;
            
            % 1.3 Fill the preallocated vector using clean block tracking
            current_idx = 3;
            for n = 1:num_neurons
                % Extract local spikes
                spikes_A = CellMatrix{idx_selected(n), st_a, rp_a};
                spikes_B = CellMatrix{idx_selected(n), st_b, rp_b};
                
                len_A = length(spikes_A);
                len_B = length(spikes_B);
                
                if len_A > 0
                    all_spikes_pair(current_idx : current_idx + len_A - 1) = spikes_A;
                    current_idx = current_idx + len_A;
                end
                if len_B > 0
                    all_spikes_pair(current_idx : current_idx + len_B - 1) = spikes_B;
                    current_idx = current_idx + len_B;
                end
            end
            
            % 1.4 Generate the unique master time profile with strict preallocation
            t_unique = unique(all_spikes_pair);
            t_filtered = t_unique([true, diff(t_unique) > 1e-9]); % Filter out numerical duplicates
            
            % Determine the exact final size and preallocate t_master
            len_filtered = length(t_filtered);
            if t_filtered(end) < tmax
                t_master = zeros(1, len_filtered + 1);
                t_master(1:len_filtered) = t_filtered;
                t_master(end) = tmax; % Append tmax without dynamic concatenation
            else
                t_master = t_filtered;
            end
            
            % --- STEP 2: Preallocate and compute the instantaneous distance profile for each labeled line ---
            % Preallocation of the population profile matrix for the current time grid size
            S_population_profiles = zeros(num_neurons, length(t_master));
            
            for n = 1:num_neurons
                train_A = sort(CellMatrix{idx_selected(n), st_a, rp_a});
                train_B = sort(CellMatrix{idx_selected(n), st_b, rp_b});
                
                if isempty(train_A) && isempty(train_B)
                    S_population_profiles(n, :) = 0;
                elseif isempty(train_A) || isempty(train_B)
                    S_population_profiles(n, :) = 1;
                else
                    % Map master time points to local spike intervals using histcounts
                    [~, ~, bin_A] = histcounts(t_master, [-Inf, train_A, Inf]);
                    [~, ~, bin_B] = histcounts(t_master, [-Inf, train_B, Inf]);
                    
                    % Local interval boundaries for Train A
                    idx_p_A = max(1, bin_A - 1); idx_n_A = min(length(train_A), bin_A);
                    x_p = train_A(idx_p_A); x_p(bin_A - 1 < 1) = tmin;
                    x_a = train_A(idx_n_A); x_a(bin_A > length(train_A)) = tmax;
                    isi_x = x_a - x_p; isi_x(isi_x <= 0) = 1e-9; % Protect against division by zero
                    
                    % Local interval boundaries for Train B
                    idx_p_B = max(1, bin_B - 1); idx_n_B = min(length(train_B), bin_B);
                    y_p = train_B(idx_p_B); y_p(bin_B - 1 < 1) = tmin;
                    y_a = train_B(idx_n_B); y_a(bin_B > length(train_B)) = tmax;
                    isi_y = y_a - y_p; isi_y(isi_y <= 0) = 1e-9; % Protect against division by zero
                    
                    % Instantaneous time differences to bounding spikes
                    dt_x_p = t_master - x_p; dt_x_a = x_a - t_master;
                    dt_y_p = t_master - y_p; dt_y_a = y_a - t_master;
                    
                    % Vectorized nearest neighbor search within the same dedicated channel
                    idx_near_B = knnsearch(train_B', t_master');
                    min_dxy = abs(t_master - train_B(idx_near_B));
                    
                    idx_near_A = knnsearch(train_A', t_master');
                    min_dyx = abs(t_master - train_A(idx_near_A));
                    
                    % Symmetric local distance contributions
                    S_1 = (min_dyx .* dt_x_a + min_dyx .* dt_x_p) ./ isi_x;
                    S_2 = (min_dxy .* dt_y_a + min_dxy .* dt_y_p) ./ isi_y;
                    
                    % Instantaneous SPIKE-distance profile for neuron n
                    denom = (isi_x + isi_y).^2;
                    S_population_profiles(n, :) = (2 * (S_1 .* isi_y + S_2 .* isi_x)) ./ denom;
                end
            end
            
            % --- STEP 3: Instantaneous population averaging (The core of Labeled Line) ---
            S_t_ll_combined = mean(S_population_profiles, 1);
            
            % --- STEP 4: Final temporal integration over the time window ---
            dval = trapz(t_master, S_t_ll_combined) / (tmax - tmin);
            
            % Save calculated distance to symmetric matrix indices
            MatrixD_LL(t_a, t_b) = dval;
            MatrixD_LL(t_b, t_a) = dval;
        end
    end
 
    %% 2. Global Decoding Performance Verification (P_LL)
    stim_A = floor(((1:num_trials)-1)/R) + 1;
    [Grid_A, Grid_B] = meshgrid(stim_A, stim_A);
    tri_upper = triu(true(num_trials), 1);  
    is_intra = (Grid_A == Grid_B) & tri_upper;
    is_inter = (Grid_A ~= Grid_B) & tri_upper;
    
    % Assign decoding performance index to the official return variable
    P_LL = mean(MatrixD_LL(is_inter)) - mean(MatrixD_LL(is_intra));
    if isnan(P_LL), P_LL = -Inf; end
end