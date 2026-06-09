%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function candidate = random_neighbor(population,N)

    candidate = population;

    nPop = length(population);

    %% ---------------------------------
    %% population complète
    %% ---------------------------------

    if nPop == N

        idx = randi(N);

        candidate(idx) = [];

        return
    end

    %% ---------------------------------
    %% population taille 1
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