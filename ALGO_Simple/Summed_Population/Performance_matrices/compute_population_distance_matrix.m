%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function D = compute_population_distance_matrix(trials,Tmax,Distances,threshold)
% COMPUTE_POPULATION_DISTANCE_MATRIX Computes the symmetric pairwise distance matrix across trials.
%
%   Calculates the pairwise spike train distance matrix :math:`D` of dimensions `[T x T]` 
%   between all pooled trials. The function automatically adds auxiliary boundary spikes 
%   to handle edge effects at `t = 0` and `t = Tmax`, computes the adaptive MRTS threshold 
%   if set to `'auto'`, and calls `SPIKE_dist_2x2` for pairwise evaluation.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Compute distance matrix with classic SPIKE distance
%      D = compute_population_distance_matrix(trials, Tmax, [1 0 0 0], 0);
%
%      % Compute distance matrix with adaptive SPIKE distance and automatic threshold
%      D = compute_population_distance_matrix(trials, Tmax, [0 0 1 0], 'auto');
%
%   :param trials: 1D cell array of length `T` containing pooled spike time vectors for each trial.
%   :type trials: cell array
%   :param Tmax: Upper temporal boundary of the analysis window.
%   :type Tmax: double
%   :param Distances: Distance metric selection mask array (e.g., `[1 0 0 0]` for SPIKE distance).
%   :type Distances: 0 or 1 array
%   :param threshold: MRTS threshold value for adaptive distances (set `0` for classic mode or `'auto'`).
%   :type threshold: double or char
%
%   :returns: **D** -- Symmetric distance matrix of dimensions `[T x T]` where `D(i,j)` represents the distance between trial `i` and trial `j`.
%   :rtype D: matrix of doubles
%
%   :Author: Maxime BELTOISE
%   :Date: May 2026
    
    T = length(trials);
    
    D = zeros(T);
    
    %% =====================================
    %% auxiliary spikes
    %% =====================================
    
    [trials_aux,aux_begin,aux_end] = add_auxiliary_spikes(trials,0,Tmax);
    
    %% =====================================
    %% pairwise distances
    %% =====================================
    
    if not (isnumeric(threshold))
        threshold = autoMRTS(trials_aux,threshold);
    end

    for i = 1:T
    
        for j = i+1:T
    
            [d,~] = SPIKE_dist_2x2(trials_aux{i}, trials_aux{j}, 0, Tmax, aux_begin(i), aux_end(i), aux_begin(j), aux_end(j), Distances, threshold);
            
            D(i,j) = d(find(Distances));
            D(j,i) = d(find(Distances));
    
        end
    end
end


