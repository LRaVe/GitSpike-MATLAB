%% Fonction to get the Time difference Matrix
% Author : Lucas RAVELOARINORO
% Date : 06/12/26

function [time_difference_matrix] = f_TD_matrix(trains,t_min,t_max)
    % Calculate the time difference matrix for a set of spike trains
    % trains: cell array where trains{i} contains spike times for train i
    % t_min: minimum time
    % t_max: maximum time
    
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