%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function accept = ...
    metropolis_acceptance( ...
        Pcandidate,...
        Pcurrent,...
        T)

    if Pcandidate > Pcurrent

        accept = true;

        return
    end

    deltaP = Pcurrent-Pcandidate;

    q = exp(-deltaP/T);

    accept = rand < q;

end