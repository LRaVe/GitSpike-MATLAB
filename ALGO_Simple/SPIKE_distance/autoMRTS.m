%% =========================================================
% Threshold (mean of the second moments of the ISIs)
% See Measures of spike train synchrony for data with multiple time scales for more informations
% =========================================================

function [MRTS] = autoMRTS(spikes)
    sum_isi_sqr = 0;
    num_isi = 0;
    for i=1:length(spikes)
        for j=1:(length(spikes{i})-1)
            sum_isi_sqr = sum_isi_sqr + (spikes{i}(j+1)-spikes{i}(j))^2;
            num_isi = num_isi + 1;
        end
    end
    MRTS = (sum_isi_sqr/num_isi)^0.5;
end