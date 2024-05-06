function S=retract(eS,extFact)
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
% Outputs: 
%   S - the row-vise array of sampled signals

S = eS(1:(size(eS,1)/extFact),1:(size(eS,2)+1-extFact));

end