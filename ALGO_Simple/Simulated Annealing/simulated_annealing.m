%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function result = simulated_annealing( ...
                    spikes,...
                    Tmax,...
                    Distances,...
                    threshold,...
                    paramsSA)

    N = size(spikes,1);

    %% =====================================
    %% cache
    %% =====================================

    cache = containers.Map;

    %% =====================================
    %% température initiale
    %% =====================================

    T0 = initialize_temperature( ...
            spikes,...
            Tmax,...
            Distances,...
            threshold,...
            paramsSA.N0);

    fprintf('Initial temperature T0 = %.4f\n',T0);

    T = T0;

    %% =====================================
    %% population initiale
    %% =====================================

    currentPop = randperm(N,randi([1 N]));

    [Pcurrent,~,~,cache] = ...
        evaluate_population_cached( ...
            spikes,...
            currentPop,...
            Tmax,...
            Distances,...
            threshold,...
            cache);

    %% =====================================
    %% meilleur état rencontré
    %% =====================================

    bestPop = currentPop;
    bestP = Pcurrent;

    %% =====================================
    %% historique
    %% =====================================

    history.P = [];
    history.bestP = [];
    history.size = [];
    history.temperature = [];

    %% =====================================
    %% nombre d'essais par plateau
    %% =====================================

    stepsPerTemp = 10*N;

    %% =====================================
    %% boucle principale
    %% =====================================

    while true

        populationChanged = false;

        for k = 1:stepsPerTemp

            candidate = ...
                random_neighbor(currentPop,N);

            [Pcandidate,~,~,cache] = ...
                evaluate_population_cached( ...
                    spikes,...
                    candidate,...
                    Tmax,...
                    Distances,...
                    threshold,...
                    cache);

            accept = ...
                metropolis_acceptance( ...
                    Pcandidate,...
                    Pcurrent,...
                    T);

            if accept

                if ~isequal(sort(candidate),...
                            sort(currentPop))

                    populationChanged = true;

                end

                currentPop = candidate;
                Pcurrent = Pcandidate;

            end

            %% meilleur rencontré

            if Pcurrent > bestP

                bestP = Pcurrent;
                bestPop = currentPop;

            end

            %% historique

            history.P(end+1) = Pcurrent;
            history.bestP(end+1) = bestP;
            history.size(end+1) = length(currentPop);
            history.temperature(end+1) = T;

        end

        %% =====================================
        %% convergence
        %% =====================================

        if ~populationChanged

            %% article :
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

        %% refroidissement

        T = T*paramsSA.coolingFactor;

    end

    %% =====================================
    %% sortie
    %% =====================================

    result.bestPopulation = sort(bestPop);

    result.bestP = bestP;

    result.history = history;

    result.uniquePopulations = cache.Count;
    
    fprintf('Iterations : %d\n',length(history.P));

    fprintf('Unique populations evaluated : %d\n', ...
            cache.Count);

end


