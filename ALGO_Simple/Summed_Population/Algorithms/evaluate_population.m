%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function [P,D,labels] = evaluate_population(spikes,neurons,Tmax,Distances,threshold)
% EVALUATE_POPULATION Computes population distance matrix and discrimination performance.
%
%   Evaluates the classification capability of a specific subset of neurons. 
%   The function concatenates spike trains for the selected neurons, computes the 
%   pairwise distance matrix :math:`D` across all trials using the requested distance 
%   metrics, and calculates the overall discrimination performance :math:`P`.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Evaluate population performance with standard distance metric
%      [P, D, labels] = evaluate_population(spikes, [1, 3, 5], Tmax, Distances, 0);
%
%      % Evaluate full population with adaptive distance threshold
%      [P, D, labels] = evaluate_population(spikes, 1:N, Tmax, Distances, 'auto');
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike timestamps.
%   :type spikes: cell or double
%   :param neurons: Vector of 1-based indices defining the subset of neurons to evaluate.
%   :type neurons: vector of integers
%   :param Tmax: Upper temporal boundary of the analysis window.
%   :type Tmax: double
%   :param Distances: Distance metric type (e.g., `'SPIKE'`, `'RI-SPIKE'`, `'SPIKE_adaptive'`, `'RI-SPIKE adaptative'`).
%   :type Distances: 0 or 1 array (e.g., `[1 0 0 0]` for the SPIKE distance)
%   :param threshold: MRTS threshold value for adaptive distances (set `0` for classic mode or `'auto'`).
%   :type threshold: double or char
%
%   :returns: **P** -- Overall discrimination performance score :math:`P` achieved by the selected subpopulation.
%   :type P: double
%   :returns: **D** -- Pairwise distance matrix :math:`D` computed between all concatenated trials.
%   :type D: matrix
%   :returns: **labels** -- Class label vector associated with each pooled trial.
%   :type labels: vector

    %% pooled trials

    [trials,labels] = build_trials(spikes,neurons);

    %% distance matrix

    D = compute_population_distance_matrix(trials,Tmax,Distances,threshold);

    %% performance

    P = compute_discrimination_performance(D,labels);

end