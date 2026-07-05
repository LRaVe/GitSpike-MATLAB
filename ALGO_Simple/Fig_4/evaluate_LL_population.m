%% ALGO Computation
% Author: Maxime BELTOISE
% Date: July 2026


function result = evaluate_LL_population(spikes,Tmax,Distances,threshold)
    
    %======================================================================
    % Evaluate Labeled Line population
    %
    % Reimplementation following:
    % Delis et al., 2015
    %
    %======================================================================
    
    [num_neurons,num_stimuli,num_repetitions] = size(spikes);
    
    T = num_stimuli*num_repetitions;
    
    %% ==============================================================
    %% Storage
    %% ==============================================================
    
    result.DistanceMatrix  = cell(1,num_neurons);
    result.Discrimination  = cell(1,num_neurons);
    result.Performance     = cell(1,num_neurons);
    result.Mn              = cell(1,num_neurons);
    
    %% ==============================================================
    %% Compute every neuron independently
    %% ==============================================================
    
    for n = 1:num_neurons
    
        %% ----------------------------------------------------------
        % Flatten trials
        %% ----------------------------------------------------------
    
        trials = cell(1,T);
    
        idx = 1;
    
        for s = 1:num_stimuli
            for r = 1:num_repetitions
    
                x = spikes{n,s,r};
    
                if isempty(x)
                    trials{idx} = [];
                else
                    trials{idx} = x(:)';
                end
    
                idx = idx+1;
    
            end
        end
    
        %% ----------------------------------------------------------
        % Pairwise distance matrix
        %% ----------------------------------------------------------
    
        D = compute_population_distance_matrix( ...
                trials,...
                Tmax,...
                Distances,...
                threshold);
    
        result.DistanceMatrix{n} = D;
    
        %% ----------------------------------------------------------
        % Paper routine
        %% ----------------------------------------------------------
    
        [~,SMatrix,RMatrix,~,Statistics] = ...
            PerformanceValue_Eero( ...
                D,...
                num_stimuli,...
                num_repetitions);
    
        result.Discrimination{n} = SMatrix;
        result.Performance{n}    = RMatrix;
    
        %% ----------------------------------------------------------
        % Equation (11)
        %% ----------------------------------------------------------
    
        Mn = zeros(num_stimuli);
    
        for s1 = 1:num_stimuli
    
            for s2 = s1:num_stimuli
    
                inter = mean(Statistics{s1,s2});
    
                intra = mean([ ...
                    Statistics{s1,s1} ...
                    Statistics{s2,s2}]);
    
                value = SMatrix(s1,s2) * (inter-intra);
    
                Mn(s1,s2) = value;
                Mn(s2,s1) = value;
    
            end
        end
    
        result.Mn{n} = Mn;
    
    end
    
    %% ==============================================================
    %% Population performance matrix (Eq.12)
    %% ==============================================================
    
    Ppopulation = zeros(num_stimuli);
    
    BestNeuron = zeros(num_stimuli);
    
    for s1 = 1:num_stimuli
    
        for s2 = s1:num_stimuli
    
            bestValue = 0;
            bestIndex = 0;
    
            for n = 1:num_neurons
    
                value = result.Mn{n}(s1,s2);
    
                if value > bestValue
    
                    bestValue = value;
                    bestIndex = n;
    
                end
    
            end
    
            Ppopulation(s1,s2) = bestValue;
            Ppopulation(s2,s1) = bestValue;
    
            BestNeuron(s1,s2) = bestIndex;
            BestNeuron(s2,s1) = bestIndex;
    
        end
    end
    
    %% ==============================================================
    %% Optimal LL population (Eq.14)
    %% ==============================================================
    
    bestPopulation = unique(BestNeuron(:));
    bestPopulation(bestPopulation==0)=[];
    
    %% ==============================================================
    %% Global LL performance (Eq.15)
    %% ==============================================================
    
    PLL = mean(Ppopulation(Ppopulation>0));
    
    %% ==============================================================
    %% Output
    %% ==============================================================
    
    result.bestPopulation = bestPopulation;
    
    result.bestP = PLL;
    
    result.populationPerformance = Ppopulation;
    
    result.bestNeuronMatrix = BestNeuron;

end






