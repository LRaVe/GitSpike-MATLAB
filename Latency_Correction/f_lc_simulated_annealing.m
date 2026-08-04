%% Function to get the shifts using simulated annealing
% Author : Lucas RAVELOARINORO
% Date : June 2026

function [all_shifts, costs] = f_lc_simulated_annealing(trains, tmin, tmax, cost_threshold)
    % This function computes the optimal shifts for each spike train using simulated annealing to minimize the overall cost. 
    %
    % It uses a pairwise RMSE cost function to evaluate the alignment of spike trains and iteratively adjusts the shifts to find a configuration that minimizes the cost. 
    %
    % Args: 
    %     trains (cell array): A cell array where each cell contains the spike times for a specific spike train.
    %     tmin (float): The minimum time to consider for spike trains.
    %     tmax (float): The maximum time to consider for spike trains.
    %     cost_threshold (float, optional): The threshold for the cost below which the optimization will stop. If not provided, a default value of 1e-3 is used.
    %
    % Returns:
    %     all_shifts (array): An array where each element represents the optimal shift for the corresponding spike train. 
    %     costs (array): An array containing the cost values at each iteration of the simulated annealing process. 


    num_trains = length(trains);
    if nargin < 4 || isempty(cost_threshold)
        cost_threshold = 1e-3;
    end
 
    % All unique pairs (i,j) with i<j  
    pairs      = nchoosek(1:num_trains, 2);
    %num_pairs  = size(pairs, 1);
 

    train_pair_idx = cell(1, num_trains);
    for t = 1:num_trains
        train_pair_idx{t} = find(pairs(:,1) == t | pairs(:,2) == t)';
    end
 
    % ── Initial pairwise RMSE, cached coincidence offsets, and fixed correspondences
    [old_rmse, sa_old_cost, pair_diff_cells] = compute_pairwise_rmse(trains, tmin, tmax, pairs);
 
    % ── SA parameters
    T = 1;           % Starting temperature
    sim_ann_temp_fact = 100000000;
    T_end = T / sim_ann_temp_fact;
    alpha = 0.9;         % Cooling factor
 
    % ── State tracking
    current_shifts = zeros(1, num_trains); 
 
    min_cost = sa_old_cost;
    min_shifts = current_shifts;          
 
    costs = sa_old_cost;
    sum_condi  = 0;
 
    % ── Main SA loop
    while T > T_end && sa_old_cost > cost_threshold
 
        iterations = 0;
        succ_iter  = 0;
 
        % Inner loop: fixed number of attempts or successes per temperature
        while iterations < 100*num_trains && succ_iter < 10*num_trains
 
            % Pick a random train and a random displacement
            t = randi(num_trains);

            displacement = randn(1) * sa_old_cost;    % Random displacement (of the order of old_cost)
            proposed_shifts = current_shifts;
            proposed_shifts(t) = proposed_shifts(t) + displacement;

            % Recompute RMSE only for pairs that involve train t  (O(n) not O(n^2))
            new_rmse = old_rmse;
            for pac = train_pair_idx{t}
                i = pairs(pac, 1);
                j = pairs(pac, 2);
                base_diffs = pair_diff_cells{pac};
                if ~isempty(base_diffs)
                    shift_delta = proposed_shifts(j) - proposed_shifts(i);
                    new_rmse(pac) = sqrt(mean((base_diffs + shift_delta) .^ 2));
                else
                    new_rmse(pac) = tmax - tmin;
                end
            end
 
            sa_new_cost  = mean(new_rmse);
            sa_delta_cost = sa_new_cost - sa_old_cost;
 
            % Metropolis acceptance criterion
            condi = (sa_delta_cost < 0) || (exp(-sa_delta_cost / T) > rand(1));
            sum_condi = sum_condi + condi;
 
            if condi
                current_shifts    = proposed_shifts;
                old_rmse          = new_rmse;
                sa_old_cost       = sa_new_cost;
                succ_iter         = succ_iter + 1;

                % Remove any global translation so train 1 stays aligned with its original timing.
                current_shifts = current_shifts - current_shifts(1);
 
                % Track the best state seen
                if sa_new_cost < min_cost
                    min_cost   = sa_new_cost;
                    min_shifts = current_shifts;
                end
            end
 
            iterations = iterations + 1;
            costs(end + 1) = sa_old_cost; %#ok<AGROW>

            if sa_old_cost <= cost_threshold
                break;
            end
        end
 
        T = T * alpha;   % Cool down
 
        if succ_iter == 0 || sa_old_cost <= cost_threshold
            break;       % No progress at this temperature → stop early
        end
    end
 
    % Normalize the best solution once more so train 1 has zero net shift, to reduce global shift.
    min_shifts = min_shifts - min_shifts(1);

    % f_plot_trains_with_correction applies: trains_corrected{i} = trains{i} - shifts(i)
    % so we negate the optimized relative shifts before returning them.
    all_shifts = -min_shifts;

end
 
 
%% Function to compute pairwise RMSE and coincidence offsets

function [rmse_vec, mean_cost, pair_diff_cells] = compute_pairwise_rmse(trains, tmin, tmax, pairs)
    % This function computes the pairwise RMSE for all unique pairs of spike trains, along with the mean cost and the coincidence offsets for each pair. 
    %
    % Args:
    %     trains (cell array): A cell array where each cell contains the spike times for a specific spike train. 
    %     tmin (float): The minimum time to consider for spike trains.
    %     tmax (float): The maximum time to consider for spike trains.
    %     pairs (matrix): A (num_pairs,2) matrix where each row contains a unique pair of spike train indices (i,j) with i<j.
    %
    % Returns:
    %     rmse_vec (array): An array where each element represents the RMSE for the corresponding pair of spike trains.
    %     mean_cost (float): The overall cost value, calculated as the mean of the RMSE values for all pairs.
    %     pair_diff_cells (cell array): A cell array where each cell contains the coincidence offsets for the corresponding pair of spike trains.


    num_pairs = size(pairs, 1);
    rmse_vec  = zeros(1, num_pairs);
    pair_diff_cells = cell(1, num_pairs);
    for pac = 1:num_pairs
        i = pairs(pac, 1);
        j = pairs(pac, 2);
        [C, spike_times, coincidence_times] = f_spike_synchro(trains{i}, trains{j}, tmin, tmax);
        if sum(C) > 0
            diffs = coincidence_times(C == 1) - spike_times(C == 1);
            pair_diff_cells{pac} = diffs(:);
            rmse_vec(pac) = sqrt(mean(diffs .^ 2));
        else
            pair_diff_cells{pac} = [];
            rmse_vec(pac) = tmax - tmin;
        end
    end
    mean_cost = mean(rmse_vec);
end
