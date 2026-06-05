function [trains] = f_synfire(tmin,tmax,n_trains,n_total_events,n_correct,n_inversed,overlap,step_in_cascade)
    % Create synfire trains with specified parameters
    % tmin: minimum time for spike generation
    % tmax: maximum time for spike generation
    % n_trains: number of spike trains to generate
    % n_total_events: total number of events (spikes) across all trains
    % n_correct: number of correctly ordered events (synchronized spikes)
    % n_inversed: number of inversed events (anti-synchronized spikes)
    % overlap: [0,1] indicating whether to allow overlapping spikes (<0.5 for no overlap)
    % step_in_cascade: time difference between spikes in a cascade (for both correct and inversed events)
    
    n_random = n_total_events - n_correct - n_inversed; % Number of random events

    trains = cell(1, n_trains); % Initialize cell array for spike trains
    for i = 1:n_trains
        trains{i} = zeros(1, n_total_events); % Preallocate spike times for each train
    end
    
    % Calculate cascade duration (time span from first to last spike in a cascade)
    cascade_duration = (n_trains - 1) * step_in_cascade;
    
    % Available time for event spacing (subtract cascade duration so last spike stays within tmax)
    available_time = (tmax - tmin) - cascade_duration;
    
    if overlap > 0
        % calculate the distance between events
        distance_between_events = available_time / (n_total_events - 1 + overlap); 
    else
        distance_between_events = available_time / (n_total_events - 1);
    end
    % Correctly ordered events (synchronized spikes)
    for event_id = 1:n_correct
        event_time = tmin + (event_id-1) * distance_between_events; % Calculate event time based on distance
        for train_id = 1:n_trains
            spike_time = event_time + (train_id-1) * step_in_cascade;
            trains{train_id}(event_id) = spike_time; % Add the synchronized spike to each train
        end
    end

    % Random events
    for event_id = 1:n_random
        event_time = tmin + n_correct * distance_between_events + (event_id-1) * distance_between_events; % Calculate event time for random events
        for train_id = 1:n_trains
            % add jitter but around the event time, and ensure the random spike time is not overlapping with the synchronized spikes
            random_spike_time = event_time + randn * 2; % Add Gaussian noise
            while any(abs(random_spike_time - trains{train_id}(1:n_correct)) < 1) % Check for overlap with synchronized spikes
                random_spike_time = event_time + randn * 2; % Regenerate if overlapping
            end
            % Clamp spike time within [tmin, tmax]
            random_spike_time = max(tmin, min(tmax, random_spike_time));
            trains{train_id}(event_id + n_correct) = random_spike_time; % Add the random spike to each train
        end
    end

    % Inversed events (anti-synchronized spikes), the first spike of the event is in the last train and the last spike is in the first
    for event_id = 1:n_inversed
        event_time = tmin + (n_correct + n_random) * distance_between_events + (event_id-1) * distance_between_events; % Calculate event time for inversed events
        for cascade_pos = 1:n_trains
            % cascade_pos=1 fires first in the last train, cascade_pos=n_trains fires last in the first train
            spike_time = event_time + (cascade_pos-1) * step_in_cascade;
            train_id = n_trains - cascade_pos + 1;  % Map cascade position to train (inverse order)
            trains{train_id}(event_id + n_correct + n_random) = spike_time; % Add the inversed spike to each train
        end
    end
    
end