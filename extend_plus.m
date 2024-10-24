function eS=extend_plus(S,extFact,S_minus,S_plus)
%EXTEND_PLUS Extends the signals (rows) in array of signals S by adding the delayed repetitions of each signal (row). 
%
% eS=extend_plus(S,extFact,S_minus,S_plus);
%
% Extends the signals (rows) in array of signals S by adding the delayed
% repetitions of each signal (row).
%
% Inputs:
%   S - the row-vise array of sampled signals
%   extFact - the number of delayed repetitions of each signal (row) in S
%   S_minus - Part that is subbed in for zeros to pad beginning of signal
%   S_plus - The part that is subbed in for zeros to pad end of signal
% Outputs: 
%   eS - extended row-vise array of signals
%
%   Adapted by Max Murphy from original code `extend` by:
% AUTHOR: Ales Holobar, FEECS, University of Maribor, Slovenia

[r,c]=size(S);
eS=zeros(r*extFact,c+extFact-1);
minus_vec = fliplr(1:extFact-1);

for k=1:r
    for m=1:extFact
        eS((k-1)*extFact+m,:)=[fliplr(S_minus(k,minus_vec(1:(m-1)))) S(k,:) S_plus(k,1:(extFact-m))];
    end
end

end