function MU_ID = generate_muap_id(MUPulses, options)
%GENERATE_MUAP_ID  Generates MUAP ID (cell array of identifier char arrays) from MUPulses OR IPTs.
%
% Syntax:
%   MU_ID = ckc.generate_muap_id(MUPulses,'Name', value, ...);
%       OR
%   MU_ID = ckc.generate_muap_id(IPTs,'Name',value,...);
%
% Inputs:
%   MUPulses (1 x nMUAPs cell) or IPTs (nMUAPs x nSamples array)
%
% Options:
%   'Prefix' {mustBeTextScalar} = "A"; % Goes before the numeric part of ID
%
% Output:
%   MU_ID - 1 x nMUAPs cell of IDs like "A-02" etc.
%
% See also: Contents    

arguments
    MUPulses % (or IPTs)
    options.Prefix {mustBeTextScalar} = "A";
end

if iscell(MUPulses)
    nMUAPs = numel(MUPulses);
else
    nMUAPs = size(MUPulses, 1); % Then we passed IPTs instead of MUPulses
end
MU_ID = cell(1, nMUAPs);
for ii = 1:nMUAPs
    MU_ID{ii} = sprintf("%s-%02d",options.Prefix,ii);
end

end