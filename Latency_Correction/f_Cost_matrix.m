%% Function to get the Cost Matrix
% Author : Lucas RAVELOARINORO
% Date : 06/12/26

function [Cost_matrix,Cost_value] = f_Cost_matrix(trains)
    % Calculate the time difference matrix for a set of spike trains
    % trains: cell array where trains{i} contains spike times for train i
    % n_events: number of events to consider
    addpath('SPIKE_synchro');
    addpath('spike_common');

    n_trains = length(trains);
    Cost_matrix = zeros(n_trains,n_trains); % Initialize numeric matrix for time differences
    
    for i=1:n_trains
        for j=[1:i-1, i+1:n_trains] % Loop through pairs of trains (excluding diagonal)
            [C, spike_times, coincidence_times] = f_spike_synchro(trains{i}, trains{j}, min(trains{i}), max(trains{i})); % Get coincidence vector and corresponding spike times
            if sum(C)~=0 % If there are coincident spikes
                for k=1:length(C)
                    if C(k) == 1 % If the spikes are coincident
                        Cost_matrix(i,j) = Cost_matrix(i,j) + (coincidence_times(k) - spike_times(k)).^2; % Calculate cost
                    end
                end
                Cost_matrix(i,j) = sqrt(Cost_matrix(i,j) / sum(C)); % Average cost for coincident spikes
            end
            if abs(Cost_matrix(i,j)) < 1e-10 % Handle numerical precision issues
                Cost_matrix(i,j) = 0;
            end
        end
    end
    % Calculate the overall cost value as the mean of the upper triangle of the Cost_matrix (excluding diagonal)
    Cost_value = mean(Cost_matrix(triu(true(size(Cost_matrix)), 1)));
end