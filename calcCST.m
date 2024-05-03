function [cst,filt_cst]=calcCST(MUPulses,fsamp,sigLen)
% INPUTS:
%  - MUPulses: cell array with MU firing times (in samples) in each cell (MUPulses{1} stores the firing times of 1st MU)
%  - fsamp: sampling frequency
%  - sigLen: length of the signals (the last MU firings does not denote the signal length).
% OUTPUTS:
%  - cst: cummulative spike train (not filtered)
%  - filt_cst: filtered cummulative spike train (hann window of length of 400 ms)
%
% PARAMETERS:
HannWinLength = 0.4; % hann window of length of 400 ms

% calculate CST
cst = zeros(1,4000+sigLen);
for muid = 1:length(MUPulses)
    cst(MUPulses{muid}+2000) = 1;
end

% filter CST
a = 1; b = hann(round(HannWinLength*fsamp));
b = b/sum(b);
%filt_cst = filtfilt(b,a,cst); 
filt_cst = filter(b,a,cst); 

filt_cst = filt_cst(2001:end-2000);

cst = cst(2001:end-2000);
cst = cst - mean(cst);