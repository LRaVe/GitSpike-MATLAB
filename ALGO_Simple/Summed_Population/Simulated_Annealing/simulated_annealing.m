%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function result = simulated_annealing(spikes,Tmax,Distances,threshold,paramsSA)

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


