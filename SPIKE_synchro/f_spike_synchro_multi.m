%% Adaptive SPIKE-synchronization computation for multiple spike trains
% Author: Lucas RAVELOARINORO
% Date: May 2026


function [C_matrix, C_global, sortedTimes, sortedValues, spike_synchro_data] = f_spike_synchro_multi(st, t_min, t_max, RMTS)
    % This function computes the SPIKE-synchronization for multiple spike trains.
    % It calculates pairwise coincidence matrices and builds a profile for each train.
    % It depends on whether RMTS (Real-Time Spike Train Synchronization) is enabled or not.
    % 
    % Args:
    %     st (cell array): A cell array where each cell contains the spike times for a specific spike train.
    %     t_min (float): The minimum time for the analysis window.
    %     t_max (float): The maximum time for the analysis window.
    %     RMTS (optional value): The arbitrary threshold for the adaptive coincidence detection.
    %
    % Returns:
    %     C_matrix (matrix): A square matrix where each element (i,j) represents the coincidence value between spike train i and spike train j.
    %     C_global (float): The global SPIKE-synchronization index, calculated as the mean of the upper triangle of C_matrix (excluding diagonal).
    %     sortedTimes (vector): A vector of sorted spike times across all trains.
    %     sortedValues (vector): A vector of sorted synchronization values across all trains.
    %     spike_synchro_data (cell array): A cell array storing [time, C_value] pairs for each train comparison.

    n_trains = length(st);
    C_matrix = zeros(n_trains, n_trains);
    spike_synchro_data = cell(n_trains, n_trains); % Store [time, C_value] pairs for each train comparison

    if nargin < 4
        RMTS = 0;
    end

    % Compute pairwise coincidence matrices and build one profile per train in the same pass.
    spike_synchro_profile = [];
    for i = 1:n_trains
        train_times = unique(sort(st{i}(st{i} >= t_min & st{i} <= t_max))); 
        train_values = zeros(length(train_times), 1);
        num_contributors = 0;

        for j = 1:n_trains  % Compare train i to all other trains
            if i == j
                C_matrix(i, j) = 1;  % Perfect synchronization with itself
                spike_synchro_data{i, j} = []; % No data for self-comparison
                continue;  % Skip self-comparison
            end

            [C_ij, times_ij] = f_adapt_spike_synchro(st{i}, st{j}, t_min, t_max, RMTS);

            if isempty(times_ij)
                spike_synchro_data{i, j} = [];
                continue;
            end

            % Update C_matrix
            C_matrix(i, j) = sum(C_ij) / length(times_ij);

            % Store paired time and C_value data for plotting
            spike_synchro_data{i, j} = [times_ij(:), C_ij(:)]; % Nx2 matrix: [time, C_value]

            % Accumulate the profile for train i across all other trains
            if ~isempty(train_times)
                train_values = train_values + C_ij(:);
                num_contributors = num_contributors + 1;
            end
        end

        if num_contributors > 0 && ~isempty(train_times)
            spike_synchro_profile = [spike_synchro_profile; [train_times(:), train_values ./ num_contributors]];
        end
    end

    % calculate global SPIKE-Synchronization index C_global as the mean of the upper triangle of C_matrix (excluding diagonal)
    C_global = mean(C_matrix(triu(true(size(C_matrix)), 1)));
    if isempty(spike_synchro_profile)
        sortedTimes = [];
        sortedValues = [];
    else
        spike_synchro_profile = sortrows(spike_synchro_profile, 1);
        sortedTimes = spike_synchro_profile(:, 1);
        sortedValues = spike_synchro_profile(:, 2);
    end
end