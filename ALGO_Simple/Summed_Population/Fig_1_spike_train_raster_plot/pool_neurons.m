%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function pooled = pool_neurons(spikes,neurons,s,r)
% POOL_NEURONS Concatenates and sorts spike timestamps across a subset of neurons.
%
%   Pools together all spike arrival times from a designated subset of neurons 
%   for a specific stimulus :math:`s` and repetition :math:`r`. The resulting spike train 
%   is merged and sorted in ascending chronological order to simulate population-level activity.
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike timestamps.
%   :type spikes: cell or double
%   :param neurons: Vector of 1-based indices defining the subset of neurons to pool.
%   :type neurons: vector of integers
%   :param s: Index of the target stimulus.
%   :type s: integer
%   :param r: Index of the target repetition/trial.
%   :type r: integer
%
%   :returns: **pooled** -- Chronologically sorted vector containing all concatenated spike timestamps from the selected neurons.
%   :rtype pooled: vector of doubles
%
%   :Author: Maxime BELTOISE
%   :Date: May 2026

    %% space allocation

    totalSpikes = 0;

    for n = neurons
        totalSpikes = totalSpikes + numel(spikes{n,s,r});
    end

    pooled = zeros(1,totalSpikes);

    %% filling

    idx = 1;

    for n = neurons

        v = spikes{n,s,r};

        l = numel(v);

        pooled(idx:idx+l-1) = v;

        idx = idx + l;

    end

    pooled = sort(pooled);

end