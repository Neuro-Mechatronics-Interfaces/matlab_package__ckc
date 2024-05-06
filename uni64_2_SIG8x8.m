function SIG = uni64_2_SIG8x8(uni, options)
%UNI64_2_SIG8x8 Convert unipolar grid (64 x nSamples) to SIG (8 x 8) cell array
%
% Example:
%   data = TMSiSAGA.Poly5.read("MyData.poly5");
%   iUni = contains({data.channels.name},'UNI');
%   uni = gradient(data.samples(iUni,:));
%   SIG = ckc.uni64_2_SIG8x8(uni); % Convert to SIG cell 8x8 grid
%
% Syntax:
%   SIG = ckc.uni64_2_SIG8x8(uni,'Name',value,...);
%
% Inputs:
%   uni - Unipolar sample grid data 64 x nSamples
%   
% Options:
%   'ChannelMap' (default: 1:64)
%
% Output:
%   SIG - 8x8 cell array each with 1 x nSamples array same as input.
%
% See also: Contents

arguments
    uni (64,:) {mustBeNumeric}
    options.ChannelMap {mustBePositive, mustBeInteger} = 1:64; % SIG is populated in column-major order; this is the index within uni used for sequential indexing to SIG
end
SIG = cell(8, 8);
for iCh = 1:64
    SIG{iCh} = uni(options.ChannelMap(iCh),:);
end
end