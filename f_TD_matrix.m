

function [time_difference_matrix] = f_TD_matrix(trains)
    % Calculate the time difference matrix for a set of spike trains
    % trains: cell array where trains{i} contains spike times for train i
    % n_events: number of events to consider
    addpath('SPIKE_synchro');
    addpath('spike_common');

    n_trains = length(trains);
    time_difference_matrix = zeros(n_trains,n_trains); % Initialize numeric matrix for time differences

    C_matrix = f_spike_synchro_multi(trains, 0, max(cellfun(@max, trains))); % Get the C_matrix for the trains
    
    for i=1:n_trains
        for j=1:n_trains
            if i ~= j
                for k=1:length(trains{i})
                    [min_dist, closest_j] = min(abs(trains{j} - trains{i}(k)));
                    if C_matrix(i,j) == 1 % If the spikes are coincident
                        time_difference_matrix(i,j) = time_difference_matrix(i,j) + min_dist; % Accumulate time differences
                    end
                end
            else
                time_difference_matrix(i,j) = 0; % Time difference with itself is zero
            end
            time_difference_matrix(i,j) = time_difference_matrix(i,j) / sum(C_matrix(i,j)); % Average time difference for coincident spikes
        end
    end

end