function S=retract(eS,extFact,offset)
%RETRACT Retracts the extended matrix back to original form. 
%
% S=ckc.retract(eS,extFact);
%
% Extends the signals (rows) in array of signals S by adding the delayed
% repetitions of each signal (row).
%
% Inputs:
%   eS - extended row-vise array of signals 
%   extFact - the number of delayed repetitions of each signal (row) in S
%   offset - Offset (default: 0) in terms of number of samples offset to retract to.
% Outputs: 
%   S - the row-vise array of sampled signals

arguments
    eS
    extFact
    offset = 0;
end
S = eS((1+offset):extFact:size(eS,1),1:(size(eS,2)+1-extFact));

end