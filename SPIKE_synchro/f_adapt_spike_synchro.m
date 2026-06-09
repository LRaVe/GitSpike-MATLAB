%% Adaptive SPIKE-synchronization computation between two spike trains
% Author: Lucas RAVELOARINORO
% Date: May 2026


function [C, spike]= f_adapt_spike_synchro(train1, train2, t_min, t_max, RMTS)
    % Adaptive SPIKE-synchronization between two spike trains
    % If RMTS is empty or not provided, use regular SPIKE synchronization
    if nargin < 5 || isempty(RMTS)
        f_spike_synchro(train1, train2, t_min, t_max);
        return;
    end
    
    train1_sliced = unique(sort(train1(train1 >= t_min & train1 <= t_max)));
    train2_sliced = unique(sort(train2(train2 >= t_min & train2 <= t_max)));

    C = zeros(length(train1_sliced), 1);
    spike = train1_sliced;  % Corresponding spike times for C
    if isempty(train1_sliced) || isempty(train2_sliced)
        return;
    end

    for i=1:length(train1_sliced)
        spike1 = train1_sliced(i);
        tau = f_interval(train1, spike1, t_min, t_max);
        [x_ip, x_if] = compute_ISI(train1, spike1, t_min, t_max);
        tau_ip = min(max(tau, 1/4*RMTS),1/2*x_ip);
        tau_if = min(max(tau, 1/4*RMTS),1/2*x_if);
        for j=1:length(train2_sliced)
            spike2 = train2_sliced(j);
            tau2 = f_interval(train2, spike2, t_min, t_max);
            [x_jp, x_jf] = compute_ISI(train2, spike2, t_min, t_max);
            tau_jp = min(max(tau2, 1/4*RMTS),1/2*x_jp);
            tau_jf = min(max(tau2, 1/4*RMTS),1/2*x_jf);
            if spike1 <= spike2
                tau_ij = min(tau_if, tau_jp);
            else
                tau_ij = min(tau_ip, tau_jf);
            end
            if abs(spike1 - spike2) < tau_ij
                C(i) = 1;
                break;
            else
                C(i) = 0;
            end
        end
    end







end


function tau = f_interval(spike_train, spike, t_min, t_max)
    % Calculate tau for adaptive coincidence detection
    % tau = min(forward_ISI, backward_ISI) / 2
    
    spike_index = find(abs(spike_train - spike) < 1e-10, 1);
    
    if isempty(spike_index)
        tau = 0;
        return;
    end
    
    % Distance to previous spike
    if spike_index > 1
        prev_dist = spike_train(spike_index) - spike_train(spike_index - 1);
    else
        % For first spike: use ISI to next spike (or time window edge)
        if spike_index < length(spike_train)
            prev_dist = spike_train(spike_index + 1) - spike_train(spike_index);
        else
            prev_dist = (t_max - t_min);  % Single spike fallback
        end
    end
    
    % Distance to next spike
    if spike_index < length(spike_train)
        next_dist = spike_train(spike_index + 1) - spike_train(spike_index);
    else
        % For last spike: use ISI to previous spike (or time window edge)
        if spike_index > 1
            next_dist = spike_train(spike_index) - spike_train(spike_index - 1);
        else
            next_dist = (t_max - t_min);  % Single spike fallback
        end
    end
    
    tau = min(prev_dist, next_dist) / 2;
end


function [x_ip,x_if] = compute_ISI(spike_train, spike, t_min, t_max)
    % Compute ISI before and after a given spike

    % Empty spike train
    if isempty(spike_train)
        x_ip = t_max - t_min;
        x_if = t_max - t_min;
        return;
    end

    tol = 1e-10;
    spike_index = find(abs(spike_train - spike) < tol, 1);

    if ~isempty(spike_index)
        % Spike exists in the train
        if spike_index > 1
            x_ip = spike - spike_train(spike_index - 1);
        else
            % At the first spike, reuse the forward ISI so boundary coincidences don't collapse to a zero window
            if length(spike_train) > 1
                x_ip = spike_train(spike_index + 1) - spike;
            else
                x_ip = t_max - t_min;
            end
        end

        if spike_index < length(spike_train)
            x_if = spike_train(spike_index + 1) - spike;
        else
            % At the last spike, reuse the backward ISI for the same reason
            if length(spike_train) > 1
                x_if = spike - spike_train(spike_index - 1);
            else
                x_if = t_max - t_min;
            end
        end
    else
        % Spike not exactly present: compute distance to nearest neighbors
        idx_before = find(spike_train < spike, 1, 'last');
        if isempty(idx_before)
            x_ip = spike - t_min;
        else
            x_ip = spike - spike_train(idx_before);
        end

        idx_after = find(spike_train > spike, 1, 'first');
        if isempty(idx_after)
            x_if = t_max - spike;
        else
            x_if = spike_train(idx_after) - spike;
        end
    end
end