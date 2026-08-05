%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function [P,D,labels,cache] = evaluate_population_cached(spikes,neurons,Tmax,Distances,threshold,cache)
% EVALUATE_POPULATION_CACHED Evaluates subpopulation performance with dictionary-based caching.
%
%   Wraps the `evaluate_population` calculation with a hash-map caching system (`containers.Map`) 
%   to eliminate redundant distance matrix and performance evaluations during combinatorial searches. 
%   Subpopulations are identified using a unique string key derived from their sorted neuron indices.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Initialize cache structure
%      cache.data = containers.Map();
%      cache.hits = 0;
%      cache.misses = 0;
%
%      % Evaluate population with caching mechanism
%      [P, D, labels, cache] = evaluate_population_cached(spikes, selection, Tmax, Distances, 0, cache);
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike timestamps.
%   :type spikes: cell or double
%   :param neurons: Vector of 1-based indices defining the subset of neurons to evaluate.
%   :type neurons: vector of integers
%   :param Tmax: Upper temporal boundary of the analysis window.
%   :type Tmax: double
%   :param Distances: Distance metric selection mask array (e.g., `[1 0 0 0]` for SPIKE distance).
%   :type Distances: 0 or 1 array
%   :param threshold: MRTS threshold value for adaptive distances (set `0` for classic mode or `'auto'`).
%   :type threshold: double or char
%   :param cache: Structure holding the hash map and execution statistics with fields:
%
%                 * **data** (*containers.Map*): Map storing evaluated subpopulation results keyed by sorted neuron strings.
%                 * **hits** (*integer*): Counter tracking the number of successful cache retrievals.
%                 * **misses** (*integer*): Counter tracking the number of new evaluations performed.
%   :type cache: struct
%
%   :returns: **P** -- Discrimination performance score :math:`P` for the selected subpopulation.
%   :type P: double
%   :returns: **D** -- Pairwise distance matrix of dimensions `[T x T]` between concatenated trials.
%   :type D: matrix of doubles
%   :returns: **labels** -- Class label vector associated with each pooled trial.
%   :type labels: vector of integers
%   :returns: **cache** -- Updated cache structure containing updated map entries and hit/miss counters.
%   :type cache: struct

    %% =====================================
    %% unique key of the population
    %% =====================================

    key = mat2str(sort(neurons));

    %% =====================================
    %% cache hit
    %% =====================================

    if isKey(cache.data,key)

        cache.hits = cache.hits + 1;

        tmp = cache.data(key);

        P = tmp.P;
        D = tmp.D;
        labels = tmp.labels;

        return

    end

    %% =====================================
    %% cache miss
    %% =====================================

    cache.misses = cache.misses + 1;

    [P,D,labels] = evaluate_population(spikes,neurons,Tmax,Distances,threshold);

    tmp.P = P;
    tmp.D = D;
    tmp.labels = labels;

    cache.data(key) = tmp;

end