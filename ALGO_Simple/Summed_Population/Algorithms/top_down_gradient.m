%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function result = top_down_gradient(spikes,Tmax,Distances,threshold)
% TOP_DOWN_GRADIENT Sequential backward elimination algorithm for optimal neuronal subpopulation selection.
%
%   Performs an iterative Top-Down (backward elimination) greedy search to identify 
%   the subpopulation of neurons that maximizes the classification performance :math:`P`. 
%   Starting from the full set of :math:`N` neurons, at each step, the algorithm evaluates all 
%   candidate subpopulations formed by removing one neuron and permanently eliminates the 
%   one whose removal yields the highest performance.
%
%   The computational complexity of this backward process is polynomial:
%
%   .. math::
%
%      \mathcal{O}\left(\frac{N(N+1)}{2}\right)
%
%   which drastically reduces the evaluation space compared to the :math:`2^N - 1` exhaustive search.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Run Top-Down selection with classic distance mode
%      result = top_down_gradient(spikes, Tmax, Distances, 0);
%
%      % Run Top-Down selection with adaptive distance metric
%      result = top_down_gradient(spikes, Tmax, Distances, 'auto');
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike timestamps.
%   :type spikes: cell or double
%   :param Tmax: Upper temporal boundary of the analysis window.
%   :type Tmax: double
%   :param Distances: Distance metric type (e.g., `'SPIKE'`, `'RI-SPIKE'`, `'SPIKE_adaptive'`, `'RI-SPIKE adaptative'`).
%   :type Distances: 0 or 1 array (e.g., `[1 0 0 0]` for the SPIKE distance)
%   :param threshold: MRTS threshold value for adaptive distances (set `0` for classic mode or `'auto'`).
%   :type threshold: double or char
%
%   :returns: **result** -- Structure containing search results and execution history with fields:
%
%             * **populations** (*cell*): List of candidate subpopulation vectors evaluated at each step.
%             * **P** (*vector*): Classification performances :math:`P` obtained at each step.
%             * **removedNeuron** (*vector*): Index of the neuron removed at each step.
%             * **candidateP** (*cell*): Individual candidate performances evaluated during each iteration.
%             * **candidatePop** (*cell*): Individual candidate populations evaluated during each iteration.
%             * **bestP** (*double*): Global optimal classification performance :math:`P_{opt}` achieved.
%             * **bestPopulation** (*vector*): Neuron indices corresponding to the optimal subpopulation.
%             * **bestStep** (*integer*): Step number where the optimal performance was reached.
%   :type result: struct
%
%   .. note::
%      The algorithm evaluates subpopulations down to a single remaining neuron to build 
%      the complete elimination trajectory.

    [N,~,~] = size(spikes);

    %% initial population

    currentPop = 1:N;

    %% search history

    result.populations = {};
    result.P = [];
    result.removedNeuron = [];
    
    % New data for the plot
    result.candidateP = {};
    result.candidatePop = {};

    step = 1;

    while length(currentPop) >= 1

        %% -----------------------------
        %% current population estimate
        %% -----------------------------

        [P_current,~,~] = evaluate_population(spikes,currentPop,Tmax,Distances,threshold);

        result.populations{step} = currentPop;
        result.P(step) = P_current;

        %% last population
        if isscalar(currentPop)     % faster than length(currentPop) == 1
            break;
        end

        %% -----------------------------
        %% test to remove each neuron
        %% -----------------------------

        bestP = -Inf;
        bestNeuron = NaN;
        bestPop = [];

        candidatePerf = nan(1,N);
        candidatePops = cell(1,N);

        for k = 1:length(currentPop)

            candidate = currentPop;
            candidate(k) = [];

            [P_candidate,~,~] = evaluate_population(spikes,candidate,Tmax,Distances,threshold);

            candidatePerf(currentPop(k)) = P_candidate;
            candidatePops{currentPop(k)} = candidate;

            if P_candidate > bestP

                bestP = P_candidate;
                bestNeuron = currentPop(k);
                bestPop = candidate;

            end
        end

        %% -----------------------------
        %% gradient descent
        %% -----------------------------

        result.removedNeuron(step) = bestNeuron;

        result.candidateP{step} = candidatePerf;
        result.candidatePop{step} = candidatePops;

        currentPop = bestPop;

        step = step + 1;

    end

    %% best population

    [result.bestP,idx] = max(result.P);

    result.bestPopulation = result.populations{idx};
    result.bestStep = idx;

end