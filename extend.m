function eS=extend(S,extFact)
%EXTEND Extends the signals (rows) in array of signals S by adding the delayed repetitions of each signal (row). 
%
% eS=extend(S,extFact);
%
% Extends the signals (rows) in array of signals S by adding the delayed
% repetitions of each signal (row).
%
% Inputs:
%   S - the row-vise array of sampled signals
%   extFact - the number of delayed repetitions of each signal (row) in S
% Outputs: 
%   eS - extended row-vise array of signals
%
% AUTHOR: Ales Holobar, FEECS, University of Maribor, Slovenia

[r,c]=size(S);
eS=zeros(r*extFact,c+extFact-1);

for k=1:r
    for m=1:extFact
        eS((k-1)*extFact+m,:)=[zeros(1,m-1) S(k,:) zeros(1,extFact-m)];
    end
end

% [r,c]=size(S);
% eS=zeros(r*extFact,c);
%
% for k=1:r
%     for m=0:extFact-1
%         if m == 0
%             eS((k-1)*extFact+m+1,:)=[zeros(1,m) S(k,:)];
%         else
%             tmp = diff(S(k,:),m);
%             eS((k-1)*extFact+m+1,:)=[zeros(1,m) tmp/max(tmp)*max(S(k,:))];
%         end
%     end
% end

end