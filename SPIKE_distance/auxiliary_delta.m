%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026

function delta = auxiliary_delta(spike, own_train, other_train, idx, aux_idx)
% AUXILIARY_DELTA Computes boundary-adjusted nearest-neighbor distance.
%
%   Calculates the smallest temporal distance between a target spike and all
%   spikes in another train, inheriting distances from neighboring inner spikes 
%   if the target spike is an auxiliary boundary spike.
%
%   Syntax:
%      delta = auxiliary_delta(spike, own_train, other_train, idx, aux_idx)
%
%   AUXILIARY DELTA MANAGEMENT
%
%   For auxiliary spikes:
%   - beginning auxiliary spike inherits delta from right neighbor
%   - ending auxiliary spike inherits delta from left neighbor
%
%   :param spike: Timestamp of the current spike.
%   :type spike: double
%   :param own_train: Vector of spike timestamps for the target train.
%   :type own_train: double
%   :param other_train: Vector of spike timestamps for the reference train.
%   :type other_train: double
%   :param idx: Index of the current spike within `own_train`.
%   :type idx: integer
%   :param aux_idx: Boolean flag indicating if the current spike is auxiliary.
%   :type aux_idx: logical | integer
%
%   :returns: **delta** -- Computed nearest-neighbor temporal distance.
%   :type: double

    % standard nearest-neighbor distance
    delta_std = min(abs(spike - other_train(:)));
    delta = delta_std;

    % auxiliary at beginning
    if (idx == 1) && aux_idx
        delta = min(abs(own_train(2) - other_train(:)));
    end

    % auxiliary at end
    if (idx == length(own_train)) && aux_idx
        delta = min(abs(own_train(end-1) - other_train(:)));
    end
end



