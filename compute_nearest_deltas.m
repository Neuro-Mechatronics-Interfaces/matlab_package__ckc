function Delta = compute_nearest_deltas(MUPulses, k)
%COMPUTE_NEAREST_DELTAS  Returns k x (nPulse - k) array of nearest sample distances for k nearest intervals.
%
% Syntax:
%   Delta = ckc.compute_nearest_deltas(MUPulses, k);
%
% Inputs:
%   MUPulses (1 x nPulses) - Sample instants vector for single MUAP
%   k                      - Number of nearest sample intervals to compute
%
% Output:
%   Delta - k x (nPulse - k) array where rows are arranged
%            from first to k-th neighbor (i.e. "1 away", "2 away", "3 away"). 
%
% See also: Contents

nPulse = numel(MUPulses);
Delta = zeros(k, nPulse-k);
iRow = 0;
for delta = 1:k
    iRow = iRow + 1;
    for iPulse = 1:(nPulse-k)
        Delta(iRow,iPulse) = abs(MUPulses(iPulse) - MUPulses(iPulse + delta));
    end
end
end