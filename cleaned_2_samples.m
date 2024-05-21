function samples = cleaned_2_samples(data, options)
%CLEANED_2_SAMPLES  Converts SIG to samples from loaded cleaned CKC data.
%
% Example:
%   [data,metadata,n] = io.load_cleaned_decomposition(14);
%   samples = ckc.cleaned_2_samples(data);
%
arguments
    data
    options.ApplyFiltering (1,1) logical = true;
    options.SpatialFilter {mustBeMember(options.SpatialFilter,{'None','MONO','SD Rows','SD Cols','Laplacian','CAR'})} = 'None';
    options.HPFCutoff (1,1) double = 100;
    options.SampleRate (1,1) double = 2000;
    options.ChannelMap = [];
    options.ChannelMapFile {mustBeTextScalar} = 'none';
    options.NumChannelsPerGrid (1,:) = 64;
    options.NumChannels (1,1) {mustBePositive, mustBeInteger} = 256;
end

nMin = inf;
nGrid = numel(data);
for ii = 1:nGrid
    nMin = min(nMin, numel(data(ii).SIG{1}));
end
nCh = options.NumChannels;
if isempty(options.ChannelMap)
    if strcmpi(options.ChannelMapFile,'none')
        channelMap = 1:options.NumChannels;
    else
        if contains(options.ChannelMapFile,'%d')
            chmapfile = sprintf(options.ChannelMapFile,nCh);
        else
            chmapfile = options.ChannelMapFile;
        end
        channelMap = getfield(load(chmapfile,'channelMap'),'channelMap');
    end
else
    channelMap = options.ChannelMap;
end
if size(channelMap,1)==1
    samples = nan(nMin,size(channelMap,2));
    useSampleGrid = false;
elseif size(channelMap,2)==1
    samples = nan(nMin,size(channelMap,1));
    useSampleGrid = false;
else
    samples = nan(nMin,size(channelMap,1),size(channelMap,2));
    useSampleGrid = true;
end
if isscalar(options.NumChannelsPerGrid)
    nChannelsPerGrid = ones(1,nGrid)*options.NumChannelsPerGrid;
elseif numel(options.NumChannelsPerGrid)~=nGrid
    error("If specified as an array, NumChannelsPerGrid element count (%d) must equal number of elements in `data` (%d).", numel(options.NumChannelsPerGrid), nGrid);
else
    nChannelsPerGrid = options.NumChannelsPerGrid;
end
nChannelsPerGrid = cumsum(nChannelsPerGrid);

if options.ApplyFiltering
    [b,a] = butter(3,options.HPFCutoff./(options.SampleRate/2),'high');
end

prevChannels = 0;
iGrid = 1;
for iCh = 1:nCh
    if channelMap(iCh) < 1
        continue;
    end
    [iRow,iCol] = ind2sub(size(channelMap),channelMap(iCh));
    if iCh > nChannelsPerGrid(iGrid)
        prevChannels = nChannelsPerGrid(iGrid);
        iGrid = iGrid + 1;
    end
    ch = iCh - prevChannels;
    if data(iGrid).discardChannelsVec(ch)==0
        [iSigRow,iSigCol] = ind2sub(size(data(iGrid).SIG),ch);
        if useSampleGrid
            samples(:,iRow,iCol) = filtfilt(b,a,data(iGrid).SIG{iSigRow,iSigCol}(1:nMin));
        else
            samples(:,iCh) = filtfilt(b,a,data(iGrid).SIG{iSigRow,iSigCol}(1:nMin));
        end
    end
end
if useSampleGrid
    switch options.SpatialFilter
        case 'SD Rows'
            samples = gradient(samples);
        case 'SD Cols'
            [~,samples] = gradient(samples);
        case 'Laplacian'
            samples = del2(samples);
    end
else
    base = 1:64;
    for iGrid = 1:nGrid
        vec = (iGrid-1)*64 + base;
        switch options.SpatialFilter
            case 'SD Rows'
                tmp = gradient(reshape(samples(:,vec),[],8,8));
                samples(:,vec) = reshape(tmp,[],64);
            case 'SD Cols'
                [~,tmp] = gradient(reshape(samples(:,vec),[],8,8));
                samples(:,vec) = reshape(tmp,[],64);
            case 'Laplacian'
                samples(:,vec) = reshape(del2(reshape(samples(:,vec),[],8,8)),[],64);
            case 'CAR'
                samples(:,vec) = samples(:,vec) - nanmean(samples(:,vec),2); %#ok<NANMEAN>
        end
    end
end

end