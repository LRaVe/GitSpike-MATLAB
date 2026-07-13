%% Function to compute and plot statistical discrimination matrix
% Date: July 2026
% Author : Laure WOLFF 

function  All_Matrices_M = calculate_plot_matrix_M(All_MatrixD, num_neurons, S, R, plotting)

    All_Matrices_M = zeros(S, S, num_neurons);
    stimulus = repelem(1:S, R);

    for n = 1:num_neurons
        MatrixD = All_MatrixD(:, :, n);

        % Convert NaNs (silent sweeps) to max distance 1 before stripping vectors
        MatrixD(isnan(MatrixD)) = 1;

        %% Discrimination matrix (M_n) via Wilcoxon criterion
        MatrixM = zeros(S, S);
        alpha_val = 0.001; % Value of paper strictness threshold in the 2018spaper

        for st1 = 1:S
            for st2 = (st1 + 1):S
                idx_st1 = find(stimulus == st1);
                idx_st2 = find(stimulus == st2);

                intra_1 = MatrixD(idx_st1, idx_st1);
                intra_2 = MatrixD(idx_st2, idx_st2);

                dist_intra_1 = intra_1(triu(true(length(idx_st1)), 1));
                dist_intra_2 = intra_2(triu(true(length(idx_st2)), 1));
                dist_inter = reshape(MatrixD(idx_st1, idx_st2), [], 1);

                h1 = 0; h2 = 0; h3 = 0;

                % Test 1: Inter vs Intra 1
                if ~isempty(dist_inter) && ~isempty(dist_intra_1)
                    [~, h1] = ranksum(dist_inter, dist_intra_1, 'alpha', alpha_val);
                end
                % Test 2: Inter vs Intra 2
                if ~isempty(dist_inter) && ~isempty(dist_intra_2)
                    [~, h2] = ranksum(dist_inter, dist_intra_2, 'alpha', alpha_val);
                end
                % Test 3: Intra 1 vs Intra 2
                if ~isempty(dist_intra_1) && ~isempty(dist_intra_2)
                    [~, h3] = ranksum(dist_intra_1, dist_intra_2, 'alpha', alpha_val);
                end

                % Paper Criterion Eq. 10: Discriminated if at least one test passes
                if (h1 == 1) || (h2 == 1) || (h3 == 1)
                    MatrixM(st1, st2) = 1;
                    MatrixM(st2, st1) = 1; 
                else
                    MatrixM(st1, st2) = 0;
                    MatrixM(st2, st1) = 0;
                end
            end
        end

        All_Matrices_M(:, :, n) = MatrixM;
    end

    %% --- PLOTTING SECTION ---
    if plotting == true
        stim_labels = cell(1, S);
        for st = 1:S
            stim_labels{st} = sprintf('S%d', st);
        end

        figure('Name', 'LL Discrimination Matrices Mn (Dynamic Colors)', 'Color', 'w');

        cols = ceil(sqrt(num_neurons * 1.25)); 
        rows = ceil(num_neurons / cols);

        hues = linspace(0.6, 1.6, num_neurons + 1); % Start around Blue (0.6), loop around
        hues(end) = []; % Remove duplicate wrap-around boundary
        hues = mod(hues, 1); % Keep values bounded between [0, 1]

        for n = 1:num_neurons
            subplot(rows, cols, n);
            imagesc(All_Matrices_M(:, :, n), [0 1]); 
            axis square;

            % Convert current unique hue selection to RGB
            current_rgb = hsv2rgb([hues(n), 0.9, 0.95]); 

            % Enforce original reference paper signatures for the first 3 neurons if preferred
            if num_neurons >= 3
                if n == 1, current_rgb = [0.0, 0.45, 1.0]; % Pure Blue
                elseif n == 2, current_rgb = [1.0, 0.0, 0.0]; % Pure Red
                elseif n == 3, current_rgb = [0.0, 0.65, 0.0]; % Pure Green
                end
            end

            % Form dynamic discrete map layout
            custom_map = [0, 0, 0; current_rgb];
            colormap(gca, custom_map);

            % % To plot valuesof the criterion on cells
            % for i = 1:S
            %     for j = 1:S
            %         text(j, i, num2str(All_Matrices_M(i, j, n)), 'HorizontalAlignment', 'center', ...
            %              'Color', 'w', 'FontWeight', 'bold', 'FontSize', 11);
            %     end
            % end
            % 
            set(gca, 'XTick', 1:S, 'YTick', 1:S, 'XTickLabel', stim_labels, 'YTickLabel', stim_labels, 'FontSize', 9);
            title(sprintf('Matrix M_{%d}', n ), 'FontSize', 11, 'FontWeight', 'bold');
        end
        shg;
    end
end


