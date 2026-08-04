%% ALGO Computation
% Author: Maxime BELTOISE
% Date: May 2026


function plot_distance_matrix(D,labels,titleStr)

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


