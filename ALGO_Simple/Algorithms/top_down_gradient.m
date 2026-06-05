%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function result = top_down_gradient(spikes,Tmax,Distances,threshold)

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