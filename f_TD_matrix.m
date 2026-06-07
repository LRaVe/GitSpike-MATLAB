function [time_difference_matrix] = f_TD_matrix(trains)
    % Calculate the time difference matrix for a set of spike trains
    % trains: cell array where trains{i} contains spike times for train i

    n_trains = length(trains);
    time_difference_matrix = zeros(n_trains, n_trains); % Initialize numeric matrix for time differences

    for i = 1:n_trains
        for j = 1:n_trains
            if i ~= j
                % Calculate pairwise time differences between spikes in train i and train j
                time_differences = [];
                for spike_i = trains{i}
                    for spike_j = trains{j}
                        time_differences(end+1) = abs(spike_j - spike_i); % Time difference from spike_i to spike_j
                    end
                end
                % Calculate the mean time difference for this pair
                if ~isempty(time_differences)
                    time_difference_matrix(i, j) = mean(time_differences);
                else
                    time_difference_matrix(i, j) = 0; % No spikes in one of the trains, store zero
                end
            else
                time_difference_matrix(i, j) = 0; % Self-comparison, store zero
            end
        end
    end
end