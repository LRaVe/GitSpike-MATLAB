%% ALGO Computation
% Author: Maxime BELTOISE
% Date: June 2026


function T0 = initialize_temperature(spikes,Tmax,Distances,threshold,N0)
% INITIALIZE_TEMPERATURE Computes the initial temperature for Simulated Annealing.
%
%   Performs a random walk over the subpopulation space to estimate the average 
%   performance variation :math:`\langle \Delta P \rangle` across neighboring state transitions. 
%   The initial temperature :math:`T_0` is calibrated to achieve an initial target acceptance 
%   rate of 95% for degradation moves:
%
%   .. math::
%
%      T_0 = -\frac{\langle \Delta P \rangle}{\ln(0.95)}
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Calibrate initial temperature over 50 random steps
%      T0 = initialize_temperature(spikes, Tmax, Distances, 0, 50);
%
%      % Calibrate initial temperature with adaptive threshold mode and 100 step
%      T0 = initialize_temperature(spikes, Tmax, Distances, 'auto', 100);
%
%   :param spikes: 3D cell array or matrix of dimensions `[num_neurons x num_stimuli x num_repetitions]` containing spike timestamps.
%   :type spikes: cell or double
%   :param Tmax: Upper temporal boundary of the analysis window.
%   :type Tmax: double
%   :param Distances: Distance metric selection mask array (e.g., `[1 0 0 0]` for SPIKE distance).
%   :type Distances: 0 or 1 array
%   :param threshold: MRTS threshold value for adaptive distances (set `0` for classic mode or `'auto'`).
%   :type threshold: double or char
%   :param N0: Number of random walk steps used to estimate mean performance variation.
%   :type N0: integer
%
%   :returns: **T0** -- Calibrated initial temperature :math:`T_0` ensuring a 95% initial acceptance rate.
%   :type T0: double

    N = size(spikes,1);

    %% population initiale

    pop = randperm(N,randi([1 N]));

    [Pprev,~,~] = evaluate_population(spikes,pop,Tmax,Distances,threshold);

    deltaP = zeros(1,N0);

    %% marche aléatoire

    for k = 1:N0

        pop = random_neighbor(pop,N);

        [Pnew,~,~] = evaluate_population(spikes,pop,Tmax,Distances,threshold);

        deltaP(k) = abs(Pnew-Pprev);

        Pprev = Pnew;

    end

    meanDelta = mean(deltaP);

    T0 = -meanDelta/log(0.95);

end