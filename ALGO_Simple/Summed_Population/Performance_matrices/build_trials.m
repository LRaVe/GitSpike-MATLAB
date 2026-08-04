%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function [trials,labels] = build_trials(spikes,neurons)
% BUILD_TRIALS Constructs pooled spike train trials and stimulus label vectors.
%
%   Iterates through all stimuli :math:`S` and repetitions :math:`R` to construct 
%   a 1D cell array of pooled trials for a specified subset of neurons. Each trial 
%   contains the merged spike times obtained via `pool_neurons`, accompanied by 
%   its corresponding stimulus class label.
%
%   The total number of trials :math:`T` is given by:
%
%   .. math::
%
%      T = S \times R
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike timestamps.
%   :type spikes: cell or double
%   :param neurons: Vector of 1-based indices defining the subset of neurons to pool.
%   :type neurons: vector of integers
%
%   :returns: **trials** -- 1D cell array of length :math:`T = S \times R` containing sorted spike time vectors for each trial.
%   :rtype trials: cell array
%   :returns: **labels** -- Class label vector of length :math:`T` associating each trial with its stimulus index :math:`s \in \{1, \dots, S\}`.
%   :rtype labels: vector of integers
%
%   :Author: Maxime BELTOISE
%   :Date: May 2026

    [~,S,R] = size(spikes);

    T = S*R;

    trials = cell(1,T);
    labels = zeros(1,T);

    idx = 0;

    for s = 1:S

        for r = 1:R

            idx = idx + 1;

            trials{idx} = pool_neurons(spikes,neurons,s,r);

            labels(idx) = s;

        end

    end

end

