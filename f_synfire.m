function [trains_synfire] = f_synfire(n_trains, n_spikes_per_train, t_min, t_max, jitter, shuffle)
    % Generate synfire trains with perfect synchronization, without overlap
    % n_trains: number of spike trains
    % n_spikes_per_train: number of spikes in each train
    % t_min, t_max: time window for spike generation
    % jitter: standard deviation of Gaussian noise added to spike times (0 for no jitter)
    % shuffle: proportion of trains to randomly shuffle (0 for no shuffling)
    trains_synfire = cell(1, n_trains);
    spike_times = linspace(t_min, t_max, n_spikes_per_train * n_trains); % Generate spike times evenly spaced between t_min and t_max
    
    for i =1:n_spikes_per_train * n_trains
        if jitter > 0
            noisy_spike = spike_times(i) + jitter * randn(); % Add Gaussian noise to the spike time
            spike_times(i) = max(t_min, min(t_max, noisy_spike)); % Ensure the noisy spike time is within bounds
        end
        trains_synfire{mod(i-1, n_trains)+1} = [trains_synfire{mod(i-1, n_trains)+1} spike_times(i)]; % Append the spike time
    end

    if shuffle > 0
        n_shuffled = round(shuffle * n_trains); % Calculate the number of trains to shuffle
        trains_to_shuffle = randperm(n_trains, n_shuffled); % Randomly select n_shuffled trains
    
    for train_idx = trains_to_shuffle
        % Shuffle the spike times within the selected train
        trains_synfire{train_idx} = trains_synfire{train_idx}(randperm(length(trains_synfire{train_idx})));
    end
end
