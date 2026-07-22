%% Function to compute and plot pairwise SPIKE-distance matrices 
% Date: July 2026
% Author: Laure WOLFF (codes from Maxime)

function All_Matrix_D = SPIKE_Distance_matrix(CellMatrix, num_neurons, S, R, tmin, tmax, plotting)
    num_trials = S * R;
    All_Matrix_D = zeros(num_trials, num_trials, num_neurons);

    % Required parameters to use Maxime's codes
    Distances = [1, 0, 0, 0]; 
    for n = 1:num_neurons
        raw_trains = cell(1, num_trials);
        counter = 1;
        for st = 1:S
            for rp = 1:R
                raw_trains{counter} = CellMatrix{n, st, rp};
                counter = counter + 1;
            end
        end

        threshold = autoMRTS(raw_trains);

        [prep_trains, aux_beg, aux_end] = add_auxiliary_spikes(raw_trains, tmin, tmax);

        MatrixD = zeros(num_trials, num_trials);

        %% Spike distance 
        for t_a = 1:num_trials
            for t_b = (t_a + 1):num_trials
                % Maxime's function
                [dist_vec, ~] = SPIKE_dist_2x2_matlab(...
                    prep_trains{t_a}, prep_trains{t_b}, ...
                    tmin, tmax, ...
                    aux_beg(t_a), aux_end(t_a), ...
                    aux_beg(t_b), aux_end(t_b), ...
                    Distances, threshold);

                MatrixD(t_a, t_b) = dist_vec(1); % dist_vec(1) is the SPIKE-distance
                MatrixD(t_b, t_a) = MatrixD(t_a, t_b);
            end
        end
        All_Matrix_D(:, :, n) = MatrixD;
    end

    %% plotting section
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

% %% Function to compute and plot pairwise SPIKE-distance matrices
% % Date: July 2026
% % Author: Laure WOLFF
% 
% function All_Matrix_D = SPIKE_Distance_matrix(CellMatrix, num_neurons, S, R, tmin, tmax, plotting)
%     num_trials = S * R;
%     All_Matrix_D = zeros(num_trials, num_trials, num_neurons);
% 
%     for n = 1:num_neurons
%         %% 1. Extraction des trains dans le BON ORDRE (S1-R1...S1-R5, S2-R1...)
%         raw_trains = cell(1, num_trials);
%         counter = 1;
%         for st = 1:S
%             for rp = 1:R
%                 raw_trains{counter} = sort(double(CellMatrix{n, st, rp}(:)'));
%                 counter = counter + 1;
%             end
%         end
% 
%         %% 2. Fenêtrage + spikes auxiliaires (Eq. A.1) pour chaque essai
%         prep_trains = cell(1, num_trials);
%         aux_begin   = zeros(1, num_trials);
%         aux_end     = zeros(1, num_trials);
% 
%         for i = 1:num_trials
%             train = raw_trains{i};
% 
%             if isempty(train)
%                 prep_trains{i} = [tmin, tmax];
%                 aux_begin(i) = 1;
%                 aux_end(i)   = 1;
%                 continue;
%             end
% 
%             ab = 0; ae = 0;
% 
%             if train(1) > tmin
%                 if length(train) >= 2
%                     aux1 = train(1) - max(train(1) - tmin, train(2) - train(1));
%                 else
%                     aux1 = tmin;
%                 end
%                 train = [aux1, train];
%                 ab = 1;
%             end
% 
%             if train(end) < tmax
%                 if length(train) >= 2
%                     aux2 = train(end) + max(tmax - train(end), train(end) - train(end-1));
%                 else
%                     aux2 = tmax;
%                 end
%                 train = [train, aux2];
%                 ae = 1;
%             end
% 
%             prep_trains{i} = train;
%             aux_begin(i) = ab;
%             aux_end(i)   = ae;
%         end
% 
%         MatrixD = zeros(num_trials, num_trials);
% 
%         %% 3. Boucle pairwise (Eq. 16-22)
%         for t_a = 1:num_trials
%             spikes1 = prep_trains{t_a};
%             aux1_begin = aux_begin(t_a);
%             aux1_end   = aux_end(t_a);
% 
%             for t_b = (t_a + 1):num_trials
%                 spikes2 = prep_trains{t_b};
%                 aux2_begin = aux_begin(t_b);
%                 aux2_end   = aux_end(t_b);
% 
%                 profile = zeros(0, 2); % [temps, valeur S]
% 
%                 %% --- Contribution des spikes du train 1 ---
%                 for idx_1 = 1:length(spikes1)
%                     if spikes2(1) > spikes1(idx_1)
%                         idx_2 = 1;
%                     elseif spikes2(end) <= spikes1(idx_1)
%                         idx_2 = length(spikes2) - 1;
%                     else
%                         idx_2 = find(spikes2 <= spikes1(idx_1), 1, 'last');
%                     end
% 
%                     ISI_dist_2 = spikes2(idx_2 + 1) - spikes2(idx_2);
% 
%                     % --- delta_tp_2 : distance NN de spikes2(idx_2) ---
%                     delta_tp_2 = min(abs(spikes2(idx_2) - spikes1));
%                     if (idx_2 == 1) && aux2_begin
%                         delta_tp_2 = min(abs(spikes2(2) - spikes1));
%                     end
%                     if (idx_2 == length(spikes2)) && aux2_begin
%                         delta_tp_2 = min(abs(spikes2(end-1) - spikes1));
%                     end
% 
%                     % --- delta_tf_2 : distance NN de spikes2(idx_2+1) ---
%                     delta_tf_2 = min(abs(spikes2(idx_2 + 1) - spikes1));
%                     if (idx_2 + 1 == 1) && aux2_end
%                         delta_tf_2 = min(abs(spikes2(2) - spikes1));
%                     end
%                     if (idx_2 + 1 == length(spikes2)) && aux2_end
%                         delta_tf_2 = min(abs(spikes2(end-1) - spikes1));
%                     end
% 
%                     xp_2 = spikes1(idx_1) - spikes2(idx_2);
%                     xf_2 = spikes2(idx_2 + 1) - spikes1(idx_1);
%                     S_2 = (delta_tp_2 * xf_2 + delta_tf_2 * xp_2) / ISI_dist_2;
% 
%                     if idx_1 > 1
%                         ISI_dist_1 = spikes1(idx_1) - spikes1(idx_1 - 1);
%                         S_1 = min(abs(spikes1(idx_1) - spikes2));
%                         if (idx_1 == length(spikes1)) && aux1_end
%                             S_1 = min(abs(spikes1(idx_1 - 1) - spikes2));
%                         end
%                         S = (S_1 * ISI_dist_2 + S_2 * ISI_dist_1) / (2 * (mean([ISI_dist_1, ISI_dist_2])^2));
%                         profile(end+1, :) = [spikes1(idx_1), S]; %#ok<AGROW>
%                     end
% 
%                     if idx_1 < length(spikes1)
%                         ISI_dist_1 = spikes1(idx_1 + 1) - spikes1(idx_1);
%                         S_1 = min(abs(spikes1(idx_1) - spikes2));
%                         if (idx_1 == 1) && aux1_begin
%                             S_1 = min(abs(spikes1(idx_1 + 1) - spikes2));
%                         end
%                         S = (S_1 * ISI_dist_2 + S_2 * ISI_dist_1) / (2 * (mean([ISI_dist_1, ISI_dist_2])^2));
%                         profile(end+1, :) = [spikes1(idx_1), S]; %#ok<AGROW>
%                     end
%                 end
% 
%                 %% --- Contribution des spikes du train 2 (symétrique) ---
%                 for idx_2 = 1:length(spikes2)
%                     if spikes1(1) > spikes2(idx_2)
%                         idx_1 = 1;
%                     elseif spikes1(end) <= spikes2(idx_2)
%                         idx_1 = length(spikes1) - 1;
%                     else
%                         idx_1 = find(spikes1 <= spikes2(idx_2), 1, 'last');
%                     end
% 
%                     ISI_dist_1 = spikes1(idx_1 + 1) - spikes1(idx_1);
% 
%                     delta_tp_1 = min(abs(spikes1(idx_1) - spikes2));
%                     if (idx_1 == 1) && aux1_begin
%                         delta_tp_1 = min(abs(spikes1(2) - spikes2));
%                     end
%                     if (idx_1 == length(spikes1)) && aux1_begin
%                         delta_tp_1 = min(abs(spikes1(end-1) - spikes2));
%                     end
% 
%                     delta_tf_1 = min(abs(spikes1(idx_1 + 1) - spikes2));
%                     if (idx_1 + 1 == 1) && aux1_end
%                         delta_tf_1 = min(abs(spikes1(2) - spikes2));
%                     end
%                     if (idx_1 + 1 == length(spikes1)) && aux1_end
%                         delta_tf_1 = min(abs(spikes1(end-1) - spikes2));
%                     end
% 
%                     xp_1 = spikes2(idx_2) - spikes1(idx_1);
%                     xf_1 = spikes1(idx_1 + 1) - spikes2(idx_2);
%                     S_1 = (delta_tp_1 * xf_1 + delta_tf_1 * xp_1) / ISI_dist_1;
% 
%                     if idx_2 > 1
%                         ISI_dist_2 = spikes2(idx_2) - spikes2(idx_2 - 1);
%                         S_2 = min(abs(spikes2(idx_2) - spikes1));
%                         if (idx_2 == length(spikes2)) && aux2_end
%                             S_2 = min(abs(spikes2(idx_2 - 1) - spikes1));
%                         end
%                         S = (S_1 * ISI_dist_2 + S_2 * ISI_dist_1) / (2 * (mean([ISI_dist_1, ISI_dist_2])^2));
%                         profile(end+1, :) = [spikes2(idx_2), S]; %#ok<AGROW>
%                     end
% 
%                     if idx_2 < length(spikes2)
%                         ISI_dist_2 = spikes2(idx_2 + 1) - spikes2(idx_2);
%                         S_2 = min(abs(spikes2(idx_2) - spikes1));
%                         if (idx_2 == 1) && aux2_begin
%                             S_2 = min(abs(spikes2(idx_2 + 1) - spikes1));
%                         end
%                         S = (S_1 * ISI_dist_2 + S_2 * ISI_dist_1) / (2 * (mean([ISI_dist_1, ISI_dist_2])^2));
%                         profile(end+1, :) = [spikes2(idx_2), S]; %#ok<AGROW>
%                     end
%                 end
% 
%                 %% --- Nettoyage : tri, clipping aux bornes avec interpolation, dédoublonnage ---
%                 profile = sortrows(profile, 1);
% 
%                 for i = 1:size(profile, 1)
%                     if profile(i, 1) < tmin
%                         idx = find(profile(:, 1) >= tmin, 1, 'first');
%                         profile(i, 2) = profile(i, 2) + ((profile(idx, 2) - profile(i, 2)) / (profile(idx, 1) - profile(i, 1))) * (tmin - profile(i, 1));
%                         profile(i, 1) = tmin;
%                     elseif profile(i, 1) > tmax
%                         idx = find(profile(:, 1) <= tmax, 1, 'last');
%                         profile(i, 2) = profile(idx, 2) + ((profile(i, 2) - profile(idx, 2)) / (profile(i, 1) - profile(idx, 1))) * (tmax - profile(idx, 1));
%                         profile(i, 1) = tmax;
%                     end
%                 end
% 
%                 [~, uidx] = unique(profile, 'rows', 'stable');
%                 profile = profile(uidx, :);
%                 profile = sortrows(profile, 1);
% 
%                 %% --- Intégration finale (Eq. 22) ---
%                 dval = trapz(profile(:,1), profile(:,2)) / (tmax - tmin);
% 
%                 MatrixD(t_a, t_b) = dval;
%                 MatrixD(t_b, t_a) = dval;
%             end
%         end
% 
%         All_Matrix_D(:, :, n) = MatrixD;
%     end
% 
%     %% 4. PLOTTING SECTION
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
%         cols = ceil(sqrt(num_neurons * 1.25));
%         rows = ceil(num_neurons / cols);
% 
%         for n = 1:num_neurons
%             subplot(rows, cols, n);
%             h = imagesc(All_Matrix_D(:, :, n));
%             colormap('jet');
%             axis square;
%             colorbar;
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
%             title(sprintf('Neuron %d', n), 'FontSize', 10, 'FontWeight', 'bold');
%             xlabel('Trials', 'FontSize', 8);
%             ylabel('Trials', 'FontSize', 8);
%         end
%         shg;
%     end
% end