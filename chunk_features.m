function features = chunk_features(Y, iFeatures, options)
%CHUNK_FEATURES Returns specific sample instants of extended samples.
%
% Syntax:
%   ipt = ckc.chunk_features(Y, iFeatures, options)
%
% Inputs:
%   Y        - Input signal matrix (size: nChannels x nTimesteps), where each row corresponds to a channel and 
%              each column represents a time sample.
%   iFeatures - Sample instants of features.
%   options  - Optional name-value arguments:
%       ChunkSize         - Size of data chunks to process at a time (default: 8192).
%       ExtensionFactor   - Number of extended samples to append at chunk boundaries (default: 20).
%       Verbose           - Boolean flag to display progress updates (default: true).
%
% Outputs:
%   features      - Extended samples at each instant (i.e. the motor unit pulse instants)
%
% See also: pinv, fprintf, extend_plus

arguments
    Y        % original signal nChannels x nTimesteps
    iFeatures        % Filter kernel: weights * inverse of covariance. 1 row per MU
    options.ChunkSize (1,1) {mustBePositive, mustBeInteger} = 2^13;
    options.ExtensionFactor (1,1) {mustBePositive,mustBeInteger} = 20;
    options.Verbose (1,1) logical = true;
end

[nChannels,nSamples] = size(Y);
nChunks = ceil(nSamples / options.ChunkSize);
sz = options.ChunkSize;
extFact = options.ExtensionFactor;
nFeatureSamples = numel(iFeatures);
features = zeros(nChannels*extFact,nFeatureSamples);
if options.Verbose
    fprintf(1,'Please wait, grabbing features...%03d%%\n',0);
end
iAssign = 0;
for iChunk = 1:nChunks
    vec = ((iChunk-1)*sz+1):min((iChunk*sz),nSamples);
    ii = iFeatures(ismember(iFeatures,vec))-vec(1)+1;
    if isempty(ii)
        continue;
    end
    iAssign = (iAssign(end)+1):(iAssign(end)+numel(ii));
    if iChunk == 1
        Y_minus = zeros(nChannels,extFact-1);
        Y_plus = Y(:,(iChunk*sz+1):(iChunk*sz+extFact-1));
    elseif iChunk == nChunks
        Y_minus = Y(:,((iChunk-1)*sz-extFact+2):((iChunk-1)*sz));
        Y_plus = zeros(nChannels,extFact-1);
    else
        Y_minus = Y(:,((iChunk-1)*sz-extFact+2):((iChunk-1)*sz));
        Y_plus = Y(:,(iChunk*sz+1):(iChunk*sz+extFact-1));
    end
    Ye = ckc.extend_plus(Y(:,vec),extFact,Y_minus,Y_plus);
    features(:,iAssign) = Ye(:,ii);
    if options.Verbose
        fprintf(1,'\b\b\b\b\b%03d%%\n', round(100*iChunk/nChunks));
    end
end

end