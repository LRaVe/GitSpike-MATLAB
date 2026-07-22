% % %% Script to compute the performance in summed population hypothesis 
% % % Date: May-June 2026
% % % Author : Laure WOLFF
% 
% %% Script to compute the performance in summed population hypothesis 
% % Date: May-July 2026
% % Author : Laure WOLFF (adapted with Maxime's codes)
% 
% function [P, MatrixD] = calculate_integrated_P_optimized(CellMatrix, selection, S, R, tmin, tmax)
%     idx_selected = find(selection == 1);
%     num_trials = S * R;
%     MatrixD = zeros(num_trials, num_trials);
% 
%     if isempty(idx_selected), P = -Inf; return; end
% 
%     % % Required parameters to use Maxime's codes
%     % Distances = [1, 0, 0, 0]; 
% 
%     %% 1. Pre-extraction and pooling of spike trains per trial (Summed Population)
%     raw_trains = cell(1, num_trials);
%     for t = 1:num_trials
%         st = floor((t-1)/R) + 1;
%         rp = mod(t-1, R) + 1;
% 
%         % Merge spike trains for the selected subpopulation
%         % pooled_spikes = [];
%         % for n_idx = 1:length(idx_selected)
%         %     n = idx_selected(n_idx);
%         %     pooled_spikes = [pooled_spikes, double(CellMatrix{n, st, rp}(:))'];
%         % end
%         % raw_trains{t} = sort(pooled_spikes);
%         pooled = [CellMatrix{idx_selected, st, rp}];
%         raw_trains{t} = sort(double(pooled(:))');
%     end
% 
%     %% 2. Train Preprocessing (Maxime's Thresholds and Auxiliary Spikes)
%     % threshold = 0;
%     [prep_trains, aux_beg, aux_end] = add_auxiliary_spikes(raw_trains, tmin, tmax);
% 
%     %% 3. Pairwise Distance Matrix Calculation
%     for t_a = 1:num_trials
%         for t_b = (t_a + 1):num_trials 
% 
%             % Call Maxime's function identically to the Labeled Line script
%             [dist_vec, ~] = SPIKE_dist_2x2_matlab(...
%                 prep_trains{t_a}, prep_trains{t_b}, ...
%                 tmin, tmax, ...
%                 aux_beg(t_a), aux_end(t_a), ...
%                 aux_beg(t_b), aux_end(t_b));
%                 %Distances, threshold);
% 
%             dval = dist_vec(1); % dist_vec(1) is the SPIKE-distance
% 
%             % Fill the symmetric matrix
%             MatrixD(t_a, t_b) = dval;
%             MatrixD(t_b, t_a) = dval;
%         end
%     end
% 
%     %% 4. Global Performance Verification (P)
%     stim_A = floor(((1:num_trials)-1)/R) + 1;
%     [Grid_A, Grid_B] = meshgrid(stim_A, stim_A);
%     tri_upper = triu(true(num_trials), 1);  
%     is_intra = (Grid_A == Grid_B) & tri_upper;
%     is_inter = (Grid_A ~= Grid_B) & tri_upper;
% 
%     P = mean(MatrixD(is_inter)) - mean(MatrixD(is_intra));
% 
%     % Safeguard against extreme mathematical noise or empty outputs
%     if isnan(P), P = -Inf; end
% end
% 
% function [SPIKE_distance_2x2, profile_mat] = SPIKE_dist_2x2_matlab(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end)
% 
%     % Cette version optimisée ne calcule QUE le SPIKE-distance classique (Distances(1)==1)
%     % Taille max possible du profil : 2 points par spike de chaque train
%     max_pts = 2 * (length(spikes1) + length(spikes2));
%     prof_t = zeros(max_pts, 1);
%     prof_S = zeros(max_pts, 1);
%     cnt = 0;
% 
%     %% =====================================================
%     % MAIN LOOP (contribution des spikes du train 1)
%     %% =====================================================
%     for idx_1 = 1 : length(spikes1)
% 
%         if spikes2(1) > spikes1(idx_1)
%             idx_2 = 1;
%         elseif spikes2(end) <= spikes1(idx_1)
%             idx_2 = length(spikes2) - 1;
%         else
%             idx_2 = find(spikes2 <= spikes1(idx_1), 1, 'last');
%         end
% 
%         ISI_dist_2 = spikes2(idx_2+1) - spikes2(idx_2);
%         delta_tp_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_begin);
%         delta_tf_2 = auxiliary_delta(spikes2(idx_2+1), spikes2, spikes1, idx_2+1, aux2_end);
%         xp_2 = spikes1(idx_1) - spikes2(idx_2);
%         xf_2 = spikes2(idx_2+1) - spikes1(idx_1);
%         S_2 = ((delta_tp_2*xf_2) + (delta_tf_2*xp_2)) / ISI_dist_2;
% 
%         if idx_1 > 1
%             ISI_dist_1 = spikes1(idx_1) - spikes1(idx_1-1);
%             S_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_end);
%             S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));
%             cnt = cnt + 1; prof_t(cnt) = spikes1(idx_1); prof_S(cnt) = S;
%         end
% 
%         if idx_1 < length(spikes1)
%             ISI_dist_1 = spikes1(idx_1+1) - spikes1(idx_1);
%             S_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_begin);
%             S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));
%             cnt = cnt + 1; prof_t(cnt) = spikes1(idx_1); prof_S(cnt) = S;
%         end
%     end
% 
%     %% =====================================================
%     % Contribution des spikes du train 2 (symétrique)
%     %% =====================================================
%     for idx_2 = 1 : length(spikes2)
% 
%         if spikes1(1) > spikes2(idx_2)
%             idx_1 = 1;
%         elseif spikes1(end) <= spikes2(idx_2)
%             idx_1 = length(spikes1) - 1;
%         else
%             idx_1 = find(spikes1 <= spikes2(idx_2), 1, 'last');
%         end
% 
%         ISI_dist_1 = spikes1(idx_1+1) - spikes1(idx_1);
%         delta_tp_1 = auxiliary_delta(spikes1(idx_1), spikes1, spikes2, idx_1, aux1_begin);
%         delta_tf_1 = auxiliary_delta(spikes1(idx_1+1), spikes1, spikes2, idx_1+1, aux1_end);
%         xp_1 = spikes2(idx_2) - spikes1(idx_1);
%         xf_1 = spikes1(idx_1+1) - spikes2(idx_2);
%         S_1 = ((delta_tp_1*xf_1) + (delta_tf_1*xp_1)) / ISI_dist_1;
% 
%         if idx_2 > 1
%             ISI_dist_2 = spikes2(idx_2) - spikes2(idx_2-1);
%             S_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_end);
%             S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));
%             cnt = cnt + 1; prof_t(cnt) = spikes2(idx_2); prof_S(cnt) = S;
%         end
% 
%         if idx_2 < length(spikes2)
%             ISI_dist_2 = spikes2(idx_2+1) - spikes2(idx_2);
%             S_2 = auxiliary_delta(spikes2(idx_2), spikes2, spikes1, idx_2, aux2_begin);
%             S = ((S_1*ISI_dist_2) + (S_2*ISI_dist_1)) / (2*(mean([ISI_dist_1 ISI_dist_2])^2));
%             cnt = cnt + 1; prof_t(cnt) = spikes2(idx_2); prof_S(cnt) = S;
%         end
%     end
% 
%     %% =====================================================
%     % Tronquer, trier, clipper aux bornes, dédoublonner
%     %% =====================================================
%     profile = [prof_t(1:cnt), prof_S(1:cnt)];
%     profile = sortrows(profile, 1);
% 
%     for i = 1:size(profile, 1)
%         if profile(i,1) < t_min
%             idx = find(profile(:,1) >= t_min, 1, 'first');
%             profile(i,2) = profile(i,2) + ((profile(idx,2) - profile(i,2)) / (profile(idx,1) - profile(i,1))) * (t_min - profile(i,1));
%             profile(i,1) = t_min;
%         elseif profile(i,1) > t_max
%             idx = find(profile(:,1) <= t_max, 1, 'last');
%             profile(i,2) = profile(idx,2) + ((profile(i,2) - profile(idx,2)) / (profile(i,1) - profile(idx,1))) * (t_max - profile(idx,1));
%             profile(i,1) = t_max;
%         end
%     end
% 
%     [~, uidx] = unique(profile, 'rows', 'stable');
%     profile = sortrows(profile(uidx, :), 1);
% 
%     profile_mat = {profile, [], [], []};
%     SPIKE_distance_2x2 = zeros(1, 4);
%     SPIKE_distance_2x2(1) = trapz(profile(:,1), profile(:,2)) / (t_max - t_min);
% end

%% Script to compute the performance in summed population hypothesis (exact & fast)
% Date: May-June-July 2026
% Author: Laure WOLFF

function [P, MatrixD] = calculate_integrated_P_optimized(CellMatrix, selection, S, R, tmin, tmax)
    idx_selected = find(selection == 1);
    num_trials = S * R;
    MatrixD = zeros(num_trials, num_trials);

    if isempty(idx_selected), P = -Inf; return; end

    %% Pre-extraction and pooling of spike trains per trial (Summed Population)
    raw_trains = cell(1, num_trials);
    for t = 1:num_trials
        st = floor((t-1)/R) + 1;
        rp = mod(t-1, R) + 1;
        pooled = [CellMatrix{idx_selected, st, rp}];
        raw_trains{t} = sort(double(pooled(:))');
    end

    %% Train preprocessing (Auxiliary Spikes via Maxime's function)
    [prep_trains, aux_beg, aux_end] = add_auxiliary_spikes(raw_trains, tmin, tmax);

    %% Pairwise distance matrix 
    for t_a = 1:num_trials
        st1 = prep_trains{t_a};
        a_b1 = aux_beg(t_a); 
        a_e1 = aux_end(t_a);
        
        for t_b = (t_a + 1):num_trials
            st2 = prep_trains{t_b};
            
            dval = spike_distance(st1, st2, tmin, tmax, ...
                a_b1, a_e1, aux_beg(t_b), aux_end(t_b));

            MatrixD(t_a, t_b) = dval;
            MatrixD(t_b, t_a) = dval;
        end
    end

    %% 4. Global performance 
    stim_A = floor(((1:num_trials)-1)/R) + 1;
    [Grid_A, Grid_B] = meshgrid(stim_A, stim_A);
    tri_upper = triu(true(num_trials), 1);
    is_intra = (Grid_A == Grid_B) & tri_upper;
    is_inter = (Grid_A ~= Grid_B) & tri_upper;

    P = mean(MatrixD(is_inter)) - mean(MatrixD(is_intra));
    if isnan(P), P = -Inf; end
end


function dval = spike_distance(spikes1, spikes2, t_min, t_max, aux1_begin, aux1_end, aux2_begin, aux2_end)
    n1 = length(spikes1);
    n2 = length(spikes2);

    max_pts = 2 * (n1 + n2);
    prof_t = zeros(max_pts, 1);
    prof_S = zeros(max_pts, 1);
    cnt = 0;

    %% Train 1 loop
    for idx_1 = 1:n1
        s1_val = spikes1(idx_1);
        
        if spikes2(1) > s1_val
            idx_2 = 1;
        elseif spikes2(end) <= s1_val
            idx_2 = n2 - 1;
        else
            idx_2 = find(spikes2 <= s1_val, 1, 'last');
        end

        ISI_dist_2 = spikes2(idx_2+1) - spikes2(idx_2);

        % Exact inlining of auxiliary_delta without function call overhead
        if aux2_begin && idx_2 == 1
            delta_tp_2 = min(abs(spikes2(2) - spikes1));
        else
            delta_tp_2 = min(abs(spikes2(idx_2) - spikes1));
        end

        if aux2_end && (idx_2 + 1) == n2
            delta_tf_2 = min(abs(spikes2(end-1) - spikes1));
        else
            delta_tf_2 = min(abs(spikes2(idx_2+1) - spikes1));
        end

        xp_2 = s1_val - spikes2(idx_2);
        xf_2 = spikes2(idx_2+1) - s1_val;

        S_2 = (delta_tp_2 * xf_2 + delta_tf_2 * xp_2) / ISI_dist_2;

        if idx_1 > 1
            ISI_dist_1 = s1_val - spikes1(idx_1-1);
            
            if aux1_end && idx_1 == n1
                S_1 = min(abs(spikes1(end-1) - spikes2));
            else
                S_1 = min(abs(s1_val - spikes2));
            end
            
            mean_isi = (ISI_dist_1 + ISI_dist_2) * 0.5;
            cnt = cnt + 1; 
            prof_t(cnt) = s1_val; 
            prof_S(cnt) = (S_1 * ISI_dist_2 + S_2 * ISI_dist_1) / (2 * (mean_isi^2));
        end

        if idx_1 < n1
            ISI_dist_1 = spikes1(idx_1+1) - s1_val;
            
            if aux1_begin && idx_1 == 1
                S_1 = min(abs(spikes1(2) - spikes2));
            else
                S_1 = min(abs(s1_val - spikes2));
            end
            
            mean_isi = (ISI_dist_1 + ISI_dist_2) * 0.5;
            cnt = cnt + 1; 
            prof_t(cnt) = s1_val; 
            prof_S(cnt) = (S_1 * ISI_dist_2 + S_2 * ISI_dist_1) / (2 * (mean_isi^2));
        end
    end

    %% Train 2 loop 
    for idx_2 = 1:n2
        s2_val = spikes2(idx_2);
        
        if spikes1(1) > s2_val
            idx_1 = 1;
        elseif spikes1(end) <= s2_val
            idx_1 = n1 - 1;
        else
            idx_1 = find(spikes1 <= s2_val, 1, 'last');
        end

        ISI_dist_1 = spikes1(idx_1+1) - spikes1(idx_1);

        if aux1_begin && idx_1 == 1
            delta_tp_1 = min(abs(spikes1(2) - spikes2));
        else
            delta_tp_1 = min(abs(spikes1(idx_1) - spikes2));
        end

        if aux1_end && (idx_1 + 1) == n1
            delta_tf_1 = min(abs(spikes1(end-1) - spikes2));
        else
            delta_tf_1 = min(abs(spikes1(idx_1+1) - spikes2));
        end

        xp_1 = s2_val - spikes1(idx_1);
        xf_1 = spikes1(idx_1+1) - s2_val;

        S_1 = (delta_tp_1 * xf_1 + delta_tf_1 * xp_1) / ISI_dist_1;

        if idx_2 > 1
            ISI_dist_2 = s2_val - spikes2(idx_2-1);
            
            if aux2_end && idx_2 == n2
                S_2 = min(abs(spikes2(end-1) - spikes1));
            else
                S_2 = min(abs(s2_val - spikes1));
            end
            
            mean_isi = (ISI_dist_1 + ISI_dist_2) * 0.5;
            cnt = cnt + 1; 
            prof_t(cnt) = s2_val; 
            prof_S(cnt) = (S_1 * ISI_dist_2 + S_2 * ISI_dist_1) / (2 * (mean_isi^2));
        end

        if idx_2 < n2
            ISI_dist_2 = spikes2(idx_2+1) - s2_val;
            
            if aux2_begin && idx_2 == 1
                S_2 = min(abs(spikes2(2) - spikes1));
            else
                S_2 = min(abs(s2_val - spikes1));
            end
            
            mean_isi = (ISI_dist_1 + ISI_dist_2) * 0.5;
            cnt = cnt + 1; 
            prof_t(cnt) = s2_val; 
            prof_S(cnt) = (S_1 * ISI_dist_2 + S_2 * ISI_dist_1) / (2 * (mean_isi^2));
        end
    end

    %% Boundary projection & exact integration 
    if cnt == 0, dval = 0; return; end

    profile = [prof_t(1:cnt), prof_S(1:cnt)];
    profile = sortrows(profile, 1);

    % Exact boundary clipping from original implementation
    for i = 1:size(profile, 1)
        if profile(i,1) < t_min
            idx = find(profile(:,1) >= t_min, 1, 'first');
            if ~isempty(idx) && profile(idx,1) ~= profile(i,1)
                profile(i,2) = profile(i,2) + ((profile(idx,2) - profile(i,2)) / (profile(idx,1) - profile(i,1))) * (t_min - profile(i,1));
            end
            profile(i,1) = t_min;
        elseif profile(i,1) > t_max
            idx = find(profile(:,1) <= t_max, 1, 'last');
            if ~isempty(idx) && profile(i,1) ~= profile(idx,1)
                profile(i,2) = profile(idx,2) + ((profile(i,2) - profile(idx,2)) / (profile(i,1) - profile(idx,1))) * (t_max - profile(idx,1));
            end
            profile(i,1) = t_max;
        end
    end

    [~, uidx] = unique(profile, 'rows', 'stable');
    profile = sortrows(profile(uidx, :), 1);

    dval = trapz(profile(:,1), profile(:,2)) / (t_max - t_min);
end