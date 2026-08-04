%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function candidate = random_neighbor(population,N)
% RANDOM_NEIGHBOR Generates a neighbor subpopulation by adding or removing a single neuron.
%
%   Selects an adjacent state in the subpopulation search space for local search 
%   or metaheuristics (e.g., Simulated Annealing). The function applies single-element 
%   mutations with boundary conditions:
%
%   * **Full population (** :math:`|S| = N` **):** Forces the removal of a randomly selected neuron.
%   * **Single-neuron population (** :math:`|S| = 1` **):** Forces the addition of an unselected neuron.
%   * **Intermediate sizes (** :math:`1 < |S| < N` **):** Randomly chooses between adding an available neuron or removing an existing one with equal probability (50/50).
%
%   :param population: Vector of 1-based indices defining the current neuron subpopulation.
%   :type population: vector of integers
%   :param N: Total number of neurons in the full dataset population.
%   :type N: integer
%
%   :returns: **candidate** -- Mutated subpopulation vector containing one less or one more neuron index, sorted in ascending order.
%   :rtype candidate: vector of integers
%
%   :Author: Maxime BELTOISE
%   :Date: June 2026

    candidate = population;

    nPop = length(population);

    %% ---------------------------------
    %% full population
    %% ---------------------------------

    if nPop == N

        idx = randi(N);

        candidate(idx) = [];

        return
    end

    %% ---------------------------------
    %% population size 1
    %% ---------------------------------

    if nPop == 1

        missing = setdiff(1:N,population);

        idx = missing(randi(length(missing)));

        candidate(end+1) = idx;

        candidate = sort(candidate);

        return
    end

    %% ---------------------------------
    %% add/remove 50-50
    %% ---------------------------------

    if rand < 0.5

        %% remove

        idx = randi(nPop);

        candidate(idx) = [];

    else

        %% add

        missing = setdiff(1:N,population);

        idx = missing(randi(length(missing)));

        candidate(end+1) = idx;

        candidate = sort(candidate);

    end

end