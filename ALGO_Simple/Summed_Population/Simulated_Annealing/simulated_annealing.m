%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function result = simulated_annealing(spikes,Tmax,Distances,threshold,paramsSA)
% SIMULATED_ANNEALING Optimizes neuronal subpopulation selection via Simulated Annealing.
%
%   Performs a global stochastic search to find the subpopulation of neurons 
%   that maximizes discrimination performance :math:`P`. Starting from a random 
%   subpopulation, the algorithm explores neighboring states using single-element mutations, 
%   accepting non-improving moves probabilistically according to the Metropolis criterion. 
%   Includes automatic temperature initialization (:math:`T_0`), dictionary-based evaluation 
%   caching, exponential cooling, and a re-annealing mechanism triggered upon stagnation.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Configure Simulated Annealing hyper-parameters
%      paramsSA.N0 = 50;                     % Number of steps for T0 estimation
%      paramsSA.steps = 20;                  % Iterations per temperature plateau
%      paramsSA.coolingFactor = 0.90;        % Exponential decay factor (alpha)
%      paramsSA.codingNeurons = 1:5;         % Ground-truth coding neurons (for history tracking)
%
%      % Run Simulated Annealing optimization
%      result = simulated_annealing(spikes, Tmax, Distances, 0, paramsSA);
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike timestamps.
%   :type spikes: cell or double
%   :param Tmax: Upper temporal boundary of the analysis window.
%   :type Tmax: double
%   :param Distances: Distance metric selection mask array (e.g., `[1 0 0 0]` for SPIKE distance).
%   :type Distances: 0 or 1 array
%   :param threshold: MRTS threshold value for adaptive distances (set `0` for classic mode or `'auto'`).
%   :type threshold: double or char
%   :param paramsSA: Structure containing algorithm control parameters with fields:
%
%                    * **N0** (*integer*): Number of random walk steps for initial temperature calibration.
%                    * **steps** (*integer*): Number of proposed neighbor evaluations per temperature level.
%                    * **coolingFactor** (*double*): Geometric decay factor :math:`\alpha \in (0, 1)` for cooling.
%                    * **codingNeurons** (*vector*): Indices of true coding neurons for composition tracking.
%   :type paramsSA: struct
%
%   :returns: **result** -- Structure containing optimization results and execution history with fields:
%
%             * **bestPopulation** (*vector*): Indices of neurons forming the optimal subpopulation.
%             * **bestP** (*double*): Maximum discrimination performance :math:`P_{\text{opt}}` achieved.
%             * **history** (*struct*): Detailed trajectories (`P`, `bestP`, `size`, `temperature`, `nCoding`, `nNonCoding`).
%             * **iterations** (*integer*): Total number of candidate moves evaluated.
%             * **acceptanceRate** (*double*): Overall proportion of accepted moves during the search.
%             * **cacheHits** (*integer*): Total number of redundant evaluations avoided by caching.
%             * **cacheMisses** (*integer*): Total number of unique populations evaluated.
%             * **hitRate** (*double*): Ratio of cache hits over total evaluations.
%   :rtype result: struct
%
%   .. note::
%      If no state change occurs across a full temperature plateau, the algorithm checks 
%      whether the current state is sub-optimal compared to `bestP`. If so, it triggers a 
%      **re-annealing** phase (restarting search from `bestPopulation` at `T0`); otherwise, 
%      it terminates.
%
%   .. note::
%      This algorithms has a MEX version
%
%   :Author: Maxime BELTOISE
%   :Date: June 2026
    N = size(spikes,1);

    %% =====================================
    %% cache
    %% =====================================

    cache.data = containers.Map;

    cache.hits = 0;
    cache.misses = 0;
    
    acceptedMoves = 0;
    totalMoves = 0;
    acceptedBetter = 0;
    acceptedWorse = 0;

    %% =====================================
    %% initial temperature
    %% =====================================

    T0 = initialize_temperature(spikes,Tmax,Distances,threshold,paramsSA.N0);

    fprintf('Initial temperature T0 = %.4f\n',T0);

    T = T0;

    %% =====================================
    %% initial population
    %% =====================================

    currentPop = randperm(N,randi([1 N]));

    [Pcurrent,~,~,cache] = evaluate_population_cached(spikes,currentPop,Tmax,Distances,threshold,cache);

    %% =====================================
    %% best population found
    %% =====================================

    bestPop = currentPop;
    bestP = Pcurrent;

    deltaP = [];

    %% =====================================
    %% search history
    %% =====================================

    history.P = [];
    history.bestP = [];
    history.size = [];
    history.temperature = [];
    history.acceptedWorse = [];
    history.nCoding = [];
    history.nNonCoding = [];

    %% =====================================
    %% number of tests per plateau
    %% =====================================

    stepsPerTemp = paramsSA.steps;

    %% =====================================
    %% main loop
    %% =====================================

    while true

        populationChanged = false;

        for k = 1:stepsPerTemp

            candidate = random_neighbor(currentPop,N);

            [Pcandidate,~,~,cache] = evaluate_population_cached(spikes,candidate,Tmax,Distances,threshold,cache);

            accept = metropolis_acceptance(Pcandidate,Pcurrent,T);

            totalMoves = totalMoves + 1;

            if accept

                acceptedMoves = acceptedMoves + 1;

                if Pcandidate > Pcurrent
            
                    acceptedBetter = acceptedBetter + 1;
            
                else
            
                    acceptedWorse = acceptedWorse + 1;
            
                end
            
                if ~isequal(sort(candidate),sort(currentPop))
                    populationChanged = true;
                end
            
                currentPop = candidate;
                Pcurrent = Pcandidate;
            
            end

            %% best found

            if Pcurrent > bestP

                bestP = Pcurrent;
                bestPop = currentPop;

            end

            %% search history

            history.P(end+1) = Pcurrent;
            history.bestP(end+1) = bestP;
            history.size(end+1) = length(currentPop);
            history.temperature(end+1) = T;
            deltaP(end+1) = abs(Pcandidate-Pcurrent);

            nCoding = numel(intersect(currentPop,paramsSA.codingNeurons));

            history.nCoding(end+1) = nCoding;
            
            history.nNonCoding(end+1) = length(currentPop) - nCoding;
        end

        %% =====================================
        %% convergence
        %% =====================================

        if ~populationChanged

            %% reannealing

            if bestP > Pcurrent

                fprintf('Reannealing...\n')

                currentPop = bestPop;
                Pcurrent = bestP;

                T = T0;

                continue

            else

                break

            end

        end

        %% cooling

        T = T*paramsSA.coolingFactor;

    end

    %% =====================================
    %% output
    %% =====================================

    result.bestPopulation = sort(bestPop);

    result.bestP = bestP;

    result.history = history;

    result.iterations = length(history.P);

    result.cacheHits = cache.hits;

    result.cacheMisses = cache.misses;
    
    result.acceptanceRate = acceptedMoves/totalMoves;

    hitRate = cache.hits/(cache.hits+cache.misses);

    result.hitRate = hitRate;

    fprintf('\n');
    fprintf('Iterations           : %d\n',length(history.P));

    fprintf('Acceptance rate      : %.1f %%\n',...
            100*acceptedMoves/totalMoves);
    
    fprintf('Accepted better      : %d\n',...
            acceptedBetter);
    
    fprintf('Accepted worse       : %d\n',...
            acceptedWorse);
    
    fprintf('Worse acceptance     : %.1f %%\n',...
            100*acceptedWorse/(acceptedBetter+acceptedWorse));
    
    fprintf('Cache hits           : %d\n',cache.hits);
    fprintf('Cache misses         : %d\n',cache.misses);
    fprintf('Unique populations   : %d\n',cache.misses);
    fprintf('Cache hit rate       : %.1f %%\n',100*hitRate);
    fprintf('Mean |ΔP| = %.6f\n',mean(deltaP));
    fprintf('\n');


end


