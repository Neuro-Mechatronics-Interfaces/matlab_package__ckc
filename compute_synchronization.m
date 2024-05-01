function [S, bins, N, R] = compute_synchronization(data, options)
%COMPUTE_SYNCHRONIZATION Calculate the synchronization between MUAP pairs.
%
% Example:
%   [data,metadata,n] = io.load_cleaned_decomposition(14);
%   [S, bins, N, R] = ckc.compute_synchronization(data);
%
% Outputs:
%   S - Cell array where rows are the same "trigger" MUAP and columns are
%           the same "target" MUAP pulse instants. Each cell element is a
%           histogram.
%   bins - The bin EDGES used in computing histograms in S.
%   N    - Total number of TRIGGER events (nRowsx1 vector)
%   R    - Pearson normalized correlation matrix between all pulse trains.
arguments
    data
    options.Bins (1,:) = -85:17:85;
    options.SyncLim (1,2) = [-0.1, 0.1];
end
nTotal = inf;
for ii = 1:numel(data)
    nTotal = min(size(data(ii).IPTs,2),nTotal);
end
nTrains = 0;
for ii = 1:numel(data)
    data(ii).IPTs = data(ii).IPTs(1:numel(data(ii).MUPulses),1:nTotal);
    nTrains = nTrains + numel(data(ii).MUPulses);
end
IPTs = vertcat(data.IPTs)';
Pulses = horzcat(data.MUPulses);

S = cell(nTrains,nTrains);
R = corrcoef(IPTs);
bins = options.Bins;
N = nan(nTrains,1);
for ii = 1:nTrains
    trig = Pulses{ii};
    N(ii) = numel(trig);
    for ik = (ii+1):nTrains
        counts = zeros(1,numel(bins)-1);
        targ = Pulses{ik};
        for iPulse = 1:numel(trig)
            counts = counts + histcounts(targ - trig(iPulse),bins);
        end
        S{ii,ik} = counts ./ N(ii);
    end
end

end