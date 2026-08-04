%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function accept = metropolis_acceptance(Pcandidate,Pcurrent,T)
% METROPOLIS_ACCEPTANCE Evaluates the Metropolis-Hastings acceptance criterion.
%
%   Determines whether a proposed candidate subpopulation state is accepted 
%   during Simulated Annealing optimization. Moves that improve classification 
%   performance (:math:`P_{\text{candidate}} > P_{\text{current}}`) are accepted 
%   deterministically. Non-improving moves are accepted probabilistically according 
%   to the Boltzmann distribution:
%
%   .. math::
%
%      q = \exp\left(-\frac{P_{\text{current}} - P_{\text{candidate}}}{T}\right)
%
%   :param Pcandidate: Classification performance score of the proposed neighbor subpopulation.
%   :type Pcandidate: double
%   :param Pcurrent: Classification performance score of the current subpopulation.
%   :type Pcurrent: double
%   :param T: Current system temperature in the annealing schedule.
%   :type T: double
%
%   :returns: **accept** -- Logical flag set to `true` if the candidate state is accepted, or `false` otherwise.
%   :rtype accept: logical
%
%   :Author: Maxime BELTOISE
%   :Date: June 2026

    if Pcandidate > Pcurrent

        accept = true;

        return
    end

    deltaP = Pcurrent-Pcandidate;

    q = exp(-deltaP/T);

    accept = rand < q;

end