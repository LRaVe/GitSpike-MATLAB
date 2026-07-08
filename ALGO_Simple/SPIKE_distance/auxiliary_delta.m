%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


%% =========================================================
% AUXILIARY DELTA MANAGEMENT
%
% For auxiliary spikes:
% - beginning auxiliary spike inherits delta from right neighbor
% - ending auxiliary spike inherits delta from left neighbor
%% =========================================================
function delta = auxiliary_delta(spike, own_train, other_train, idx, aux_idx)

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



