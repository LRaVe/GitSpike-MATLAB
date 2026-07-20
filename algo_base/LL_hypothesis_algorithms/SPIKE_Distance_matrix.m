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
%         Precomputed_Trains = cell(1, num_trials);
%         counter = 1;
%         for st = 1:S
%             for rp = 1:R
%                 spikes = double(CellMatrix{n, st, rp});
%                 spikes = sort(spikes(:)');
%                 spikes = spikes(spikes >= tmin & spikes <= tmax);
%                 Precomputed_Trains{counter} = spikes;
%                 counter = counter + 1;
%             end
%         end
% 
%         MatrixD = zeros(num_trials, num_trials);
% 
%         %% 2. Boucle pairwise
%         for t_a = 1:num_trials
%             train_A = Precomputed_Trains{t_a};
% 
%             for t_b = (t_a + 1):num_trials
%                 train_B = Precomputed_Trains{t_b};
% 
%                 % --- Cas limites ---
%                 if isempty(train_A) && isempty(train_B)
%                     dval = 0;
%                 elseif isempty(train_A) || isempty(train_B)
%                     dval = 1;
%                 else
%                     %% --- Construction du train A étendu (spikes auxiliaires, Eq. A.1/A.3) ---
%                     Ma = length(train_A);
%                     if Ma == 1
%                         tA_saux = tmin; tA_eaux = tmax;
%                     else
%                         tA_saux = train_A(1)   - max(train_A(1) - tmin, train_A(2) - train_A(1));
%                         tA_eaux = train_A(end) + max(tmax - train_A(end), train_A(end) - train_A(end-1));
%                     end
%                     dA_real = zeros(1, Ma);
%                     for i = 1:Ma
%                         dA_real(i) = min(abs(train_A(i) - train_B));
%                     end
%                     extA_t = [tA_saux, train_A, tA_eaux];
%                     extA_d = [dA_real(1), dA_real, dA_real(end)];
% 
%                     %% --- Construction du train B étendu ---
%                     Mb = length(train_B);
%                     if Mb == 1
%                         tB_saux = tmin; tB_eaux = tmax;
%                     else
%                         tB_saux = train_B(1)   - max(train_B(1) - tmin, train_B(2) - train_B(1));
%                         tB_eaux = train_B(end) + max(tmax - train_B(end), train_B(end) - train_B(end-1));
%                     end
%                     dB_real = zeros(1, Mb);
%                     for i = 1:Mb
%                         dB_real(i) = min(abs(train_B(i) - train_A));
%                     end
%                     extB_t = [tB_saux, train_B, tB_eaux];
%                     extB_d = [dB_real(1), dB_real, dB_real(end)];
% 
%                     %% --- Points d'évaluation et profils S_x(t), S_y(t) (Eq. 19) ---
%                     t_all = unique([tmin, train_A, train_B, tmax]);
%                     Nt = length(t_all);
% 
%                     Sx = zeros(1, Nt);
%                     Sy = zeros(1, Nt);
%                     isix = zeros(1, Nt);
%                     isiy = zeros(1, Nt);
% 
%                     for i = 1:Nt
%                         t = t_all(i);
% 
%                         j = find(extA_t <= t, 1, 'last');
%                         if j >= length(extA_t), j = length(extA_t) - 1; end
%                         xp = extA_t(j);   xa = extA_t(j + 1);
%                         dp = extA_d(j);   da = extA_d(j + 1);
%                         isix(i) = xa - xp;
%                         Sx(i) = ((t - xp) * da + (xa - t) * dp) / isix(i);
% 
%                         k = find(extB_t <= t, 1, 'last');
%                         if k >= length(extB_t), k = length(extB_t) - 1; end
%                         yp = extB_t(k);   ya = extB_t(k + 1);
%                         dpB = extB_d(k);  daB = extB_d(k + 1);
%                         isiy(i) = ya - yp;
%                         Sy(i) = ((t - yp) * daB + (ya - t) * dpB) / isiy(i);
%                     end
% 
%                     %% --- Combinaison bivariée (Eq. 20) et intégration (Eq. 22) ---
%                     denom = 0.5 * (isix + isiy) .^ 2;
%                     S_t = (Sx .* isiy + Sy .* isix) ./ denom;
%                     dval = trapz(t_all, S_t) / (tmax - tmin);
%                 end
% 
%                 MatrixD(t_a, t_b) = dval;
%                 MatrixD(t_b, t_a) = dval;
%             end
%         end
% 
%         All_Matrix_D(:, :, n) = MatrixD;
%     end
% 
%     %% 3. PLOTTING SECTION
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
