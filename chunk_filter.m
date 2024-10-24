function ipt = chunk_filter(Y, K, options)
%CHUNK_FILTER Filters the input signal using a specified kernel and returns the motor unit impulse train.
%
% Syntax:
%   ipt = ckc.chunk_filter(Y, K, options)
%
% Description:
%   This function applies a filter kernel to the input signal `Y` in chunks, producing an output impulse train
%   for each motor unit. The input signal is processed in manageable chunks to reduce memory load, with options 
%   for extending the signal and adapting the covariance matrix used in the filtering operation.
%
% Inputs:
%   Y        - Input signal matrix (size: nChannels x nTimesteps), where each row corresponds to a channel and 
%              each column represents a time sample.
%   K        - Filter kernel (size: nMUs x nChannels), where each row represents the weights for a motor unit.
%              This can either be the weights alone or weights multiplied by the inverse of the covariance matrix.
%   options  - Optional name-value arguments:
%       AdaptCovariance   - Boolean flag indicating whether to adapt the covariance matrix during filtering 
%                           (default: false).
%       ChunkSize         - Size of data chunks to process at a time (default: 8192).
%       ExtensionFactor   - Number of extended samples to append at chunk boundaries (default: 20).
%       Verbose           - Boolean flag to display progress updates (default: true).
%
% Outputs:
%   ipt      - Motor unit impulse train (size: nMUs x nTimesteps), with one row per motor unit.
%
% Notes:
%   - The input signal `Y` is processed in chunks to reduce memory usage.
%   - The signal is extended by `options.ExtensionFactor` samples at chunk boundaries to prevent edge effects.
%   - If `options.AdaptCovariance` is true, the filter kernel `K` is dynamically adapted based on the chunked data.
%     Otherwise, `K` is assumed to already include the covariance matrix.
%   - The first and last extended regions of the output impulse train are set to zero to avoid boundary artifacts.
%
% Example:
%   Y = randn(64, 10000);  % Example signal (64 channels, 10000 time points)
%   K = randn(10, 64);     % Example filter kernel for 10 motor units
%   ipt = ckc.chunk_filter(Y, K, 'AdaptCovariance', false);
%
% See also: pinv, fprintf, extend_plus

arguments
    Y        % original signal nChannels x nTimesteps
    K        % Filter kernel: weights * inverse of covariance. 1 row per MU
    options.AdaptCovariance (1,1) logical = false;
    options.ChunkSize (1,1) {mustBePositive, mustBeInteger} = 2^13;
    options.ExtensionFactor (1,1) {mustBePositive,mustBeInteger} = 20;
    options.Verbose (1,1) logical = true;
end

[nChannels,nSamples] = size(Y);
nChunks = ceil(nSamples / options.ChunkSize);
sz = options.ChunkSize;
extFact = options.ExtensionFactor;
ipt = zeros(size(K,1),size(Y,2));
if options.Verbose
    fprintf(1,'Please wait, transferring filter...%03d%%\n',0);
end
for iChunk = 1:nChunks
    vec = ((iChunk-1)*sz+1):min((iChunk*sz),nSamples);
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
    if options.AdaptCovariance % Assumes "K" is "weights" (MUAP extended waveform templates)
        P = pinv((1/sqrt(sz-1)) * (Ye * Ye')); 
        tmp = K * P * Ye;
    else % Otherwise, assume K is "weights" * P
        tmp = K * Ye;
    end
    ipt(:,vec) = tmp(:,1:numel(vec));
    if options.Verbose
        fprintf(1,'\b\b\b\b\b%03d%%\n', round(100*iChunk/nChunks));
    end
end
ipt(:,[1:(extFact+1), ((end-extFact):end)]) = 0;

end