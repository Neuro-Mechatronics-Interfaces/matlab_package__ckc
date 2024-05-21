function [snips, tsnip, pulses, samples] = cleaned_2_template(data, iGrid, iMUAP, options)
%CLEANED_2_TEMPLATE  Converts SIG to MUAP snippets template from loaded cleaned CKC data.
%
% Example:
%   [data,metadata,n] = io.load_cleaned_decomposition(14);
%   [snips, tsnip, samples] = ckc.cleaned_2_template(data); 
arguments
    data
    iGrid (1,1) {mustBePositive, mustBeInteger}
    iMUAP (1,1) {mustBePositive, mustBeInteger}
    options.ChannelMap = [];
    options.ChannelMapFile {mustBeTextScalar} = 'none';
    options.NumChannels (1,1) {mustBePositive, mustBeInteger} = 256;
    options.NumChannelsPerGrid (1,:) {mustBePositive, mustBeInteger} = 64;
    options.Vector (1,:) = -18:18;
    options.SampleRate (1,1) = 2000;
    options.LoadMethod {mustBeMember(options.LoadMethod, {'Cleaned2Samples','Sig2Samples'})} = 'Sig2Samples';
    options.SpatialFilter {mustBeMember(options.SpatialFilter,{'None','MONO','SD Rows','SD Cols','Laplacian','CAR'})} = 'None';
end
if iGrid > numel(data)
    error("iGrid (%d) must be <= number of grids (max: %d)", iGrid, numel(data));
end
if iMUAP > numel(data(iGrid).MUPulses)
    error("iMUAP (%d) must be <= number of pulse trains on grid %d (max: %d)", iMUAP, iGrid, numel(data(iGrid).MUPulses));
end
pulses = data(iGrid).MUPulses{iMUAP};
nPulse = numel(pulses);
pulses = reshape(pulses, nPulse, 1);
switch options.LoadMethod
    case 'Sig2Samples'
        samples = ckc.sig_2_samples(data);
    case 'Cleaned2Samples'
        samples = ckc.cleaned_2_samples(data,...
            'ChannelMap',options.ChannelMap, ...
            'ChannelMapFile',options.ChannelMapFile, ...
            'NumChannels',options.NumChannels, ...
            'NumChannelsPerGrid', options.NumChannelsPerGrid, ...
            'SpatialFilter', options.SpatialFilter)';
end
mask = (pulses + options.Vector)';
iExc = any((mask < 1) | (mask > size(samples,2)),1);
mask(iExc,:) = [];
nPulse = size(mask,2);
pulses(iExc) = [];
snips = nan(numel(options.Vector),size(samples,1),nPulse);
tsnip = options.Vector./(options.SampleRate*1e-3); % Return times in milliseconds

for iCh = 1:size(samples,1)
    tmp = samples(iCh,:);
    snips(:,iCh,:) = tmp(mask);
end

end