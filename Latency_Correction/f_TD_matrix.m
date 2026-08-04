%% Function to compute the time difference matrix for a set of spike trains
% Author : Lucas RAVELOARINORO
% Date : 06/12/26

function [time_difference_matrix] = f_TD_matrix(trains,t_min,t_max)
    % Calculate the time difference matrix for a set of spike trains
    % This function computes the time difference matrix for a set of spike trains. The time difference matrix is a square matrix where the element at (i,j) represents the average time difference between coincident spikes in train i and train j. 
    % 
    % Args:
    %     trains (cell array): A cell array where each cell contains the spike times for a specific spike train. 
    %     t_min (float): The minimum time to consider for spike trains. If not provided, it will be calculated as the minimum spike time across all trains. 
    %     t_max (float): The maximum time to consider for spike trains. If not provided, it will be calculated as the maximum spike time across all trains.
    %
    % Returns:
    %     time_difference_matrix (matrix): A square matrix where the element at (i,j) represents the average time difference between coincident spikes in train i and train j. If there are no coincident spikes between two trains, the corresponding element will be set to zero.

    
    if nargin < 2 || isempty(t_min)
        all_spikes = [trains{:}];
        t_min = min(all_spikes);
        t_max = max(all_spikes);
    end

    n_trains = length(trains);
    time_difference_matrix = zeros(n_trains,n_trains); % Initialize numeric matrix for time differences
    
    for i=1:n_trains
        for j=[1:i-1, i+1:n_trains] % Loop through pairs of trains (excluding diagonal)
            [C, spike_times, coincidence_times] = f_spike_synchro(trains{i}, trains{j}, t_min, t_max); % Get coincidence vector and corresponding spike times
            if sum(C)~=0 % If there are coincident spikes 
                for k=1:length(C)
                    if C(k) == 1 % If the spikes are coincident
                        time_difference_matrix(i,j) = time_difference_matrix(i,j) + (coincidence_times(k) - spike_times(k)); % Calculate time difference
                    end
                end
                time_difference_matrix(i,j) = time_difference_matrix(i,j) / sum(C); % Average time difference for coincident spikes
                if abs(time_difference_matrix(i,j)) < 1e-10 % Handle numerical precision issues
                    time_difference_matrix(i,j) = 0;
                end
            end
        end
    end
end