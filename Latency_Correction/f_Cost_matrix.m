%% Function to get the Cost Matrix
% Author : Lucas RAVELOARINORO
% Date : 06/12/26

function [Cost_matrix,Cost_value] = f_Cost_matrix(trains,t_min,t_max)
    % This function computes the cost matrix based on the time differences between coincident spikes in different spike trains.
    %
    % Args:
    %     trains (cell array): A cell array where each cell contains the spike times for a specific spike train.
    %     t_min (float): The minimum time to consider for spike trains.
    %     t_max (float): The maximum time to consider for spike trains.
    %
    % Returns:
    %     Cost_matrix (matrix): A square matrix where the element at (i,j) is the cost associated with the time differences between coincident spikes in train i and train j. 
    %     Cost_value (float): The overall cost value, calculated as the mean of the upper triangle of the Cost_matrix.

    
    if nargin < 2 || isempty(t_min)
        all_spikes = [trains{:}];
        t_min = min(all_spikes);
        t_max = max(all_spikes);
    end

    n_trains = length(trains);
    Cost_matrix = zeros(n_trains,n_trains); % Initialize numeric matrix for time differences
    
    for i=1:n_trains
        for j=[1:i-1, i+1:n_trains] % Loop through pairs of trains (excluding diagonal)
            [C, spike_times, coincidence_times] = f_spike_synchro(trains{i}, trains{j}, t_min, t_max); % Get coincidence vector and corresponding spike times
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