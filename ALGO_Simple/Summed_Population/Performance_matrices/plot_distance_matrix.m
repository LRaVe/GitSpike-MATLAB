%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function plot_distance_matrix(D,labels,titleStr)
% PLOT_DISTANCE_MATRIX Displays the pairwise trial distance matrix with stimulus separation lines.
%
%   Renders a color-mapped 2D matrix representing pairwise distances :math:`D` between trials.
%   The function automatically computes and displays the overall discrimination performance 
%   :math:`P` in the figure title, draws solid grid lines separating different stimulus classes, 
%   and dynamically adapts axis tick labels depending on the number of trials :math:`T`.
%
%   Valid call structures:
%
%   .. code-block:: matlab
%
%      % Display distance matrix with custom title
%      plot_distance_matrix(D, labels, 'SPIKE Distance Matrix');
%
%   :param D: Pairwise distance matrix of dimensions `[T x T]` computed between all trials.
%   :type D: matrix of doubles
%   :param labels: Class label vector of length `T` associating each trial with its stimulus index.
%   :type labels: vector of integers
%   :param titleStr: Title text string to display at the top of the plot.
%   :type titleStr: char or string
%
%   .. note::
%      - If the number of trials :math:`T \le 10`, all individual trial labels (`S1R1`, `S1R2`, ...) are shown.
%      - If :math:`T > 10`, an adaptive sub-sampling displays up to 10 evenly spaced tick marks to prevent label overlap.


    %% =====================================
    %% discrimination performance
    %% =====================================

    P = compute_discrimination_performance(D,labels);

    %% =====================================
    %% display matrix
    %% =====================================

    imagesc(D);

    axis square;

    colormap jet;

    colorbar;

    title(sprintf('%s   |   P = %.3f', titleStr,P), 'FontWeight','bold');

    xlabel('Recording');
    ylabel('Recording');

    hold on;
    box on;

    %% =====================================
    %% trial labels
    %% =====================================

    S = max(labels);
    R = sum(labels==1);

    T = length(labels);

    trialLabels = cell(1,T);

    idx = 0;

    for s = 1:S
        for r = 1:R

            idx = idx + 1;

            trialLabels{idx} = sprintf('S%dR%d',s,r);

        end
    end

    %% =====================================
    %% adaptive tick display
    %% =====================================

    maxLabels = 10;

    if T <= maxLabels

        tickPos = 1:T;
        tickLabels = trialLabels;

    else

        tickPos = round(linspace(1,T,maxLabels));

        % avoid duplicates caused by rounding
        tickPos = unique(tickPos);

        tickLabels = trialLabels(tickPos);

    end

    xticks(tickPos);
    yticks(tickPos);

    xticklabels(tickLabels);
    yticklabels(tickLabels);

    xtickangle(90);

    %% =====================================
    %% remove tick marks
    %% =====================================

    set(gca,'TickLength',[0 0]);

    %% =====================================
    %% separation lines of stimuli
    %% =====================================

    for s = 1:S-1

        pos = s*R + 0.5;

        xline(pos,'k','LineWidth',2);
        yline(pos,'k','LineWidth',2);

    end

end


