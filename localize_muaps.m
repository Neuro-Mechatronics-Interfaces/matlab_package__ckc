function [i_source, pk2pk, snips, tsnip] = localize_muaps(data_sync, data_cleaned, options)
%LOCALIZE_MUAPS  Returns the index to the peak "source channel" for each IPT
%
% Syntax:
%   [i_source,pk2pk,snips,tsnip] = ckc.localize_muaps(data_sync, data_cleaned, 'Name', value, ...);
%
% Output:
%   i_source - Channel index where template has largest peak-to-peak. If
%               insufficient number of MUAPs, then the corresponding index
%               value in i_source is NaN. For example if there are 7 MUAPs,
%               then i_source should be 7x1 and if MUAP-3 only has 2 spikes
%               then i_source(3) is NaN.
arguments
    data_sync    struct % Data as loaded from the synchronized single-file from TMSi
    data_cleaned struct % Data as loaded from the CKC decomposition output files
    options.ChannelOrder (1,:) = nan;
    options.MinNumberMUAPs (1,1) = 20;
    options.NumChannels (1,1) {mustBePositive, mustBeInteger} = 256;
    options.Vector (1,:) = -18:18;
    options.SampleRate (1,1) = 2000;
    options.SpatialFilter {mustBeMember(options.SpatialFilter,{'None','MONO','SD Rows','SD Cols','Laplacian','CAR'})} = 'None';
    options.Verbose (1,1) logical = true;
end
if isnan(options.ChannelOrder(1))
    channelOrder = 1:size(data_sync.uni,1);
else
    channelOrder = options.ChannelOrder;
end
N = numel(data_cleaned.MUPulses);
i_source = nan(N,1);
snips = cell(N,1);
pk2pk = cell(N,1);
if options.Verbose
    fprintf(1,'Spatially localizing MUAP template peaks...000%%\n');
end
tsnip = options.Vector ./ options.SampleRate;
nT = numel(options.Vector);
vec = options.Vector';
for ii = 1:N
    nPulse = numel(data_cleaned.MUPulses{ii});
    if nPulse < options.MinNumberMUAPs
        continue;
    end
    mask = data_cleaned.MUPulses{ii} + vec;
    snips{ii} = nan(nT, options.NumChannels, nPulse);
    for iCh = 1:options.NumChannels
        tmp = data_sync.uni(channelOrder(iCh),:);
        snips{ii}(:,iCh,:) = tmp(mask);
    end

    pk2pk{ii} = nan(options.NumChannels,1);
    for iCh = 1:options.NumChannels
        template = mean(snips{ii},3);
        pk2pk{ii}(iCh) = max(template(:,iCh)) - min(template(:,iCh));
    end

    [~,i_source(ii)] = max(pk2pk{ii});
    if options.Verbose
        fprintf(1,'\b\b\b\b\b%03d%%\n', round(100*ii/N));
    end
end
end