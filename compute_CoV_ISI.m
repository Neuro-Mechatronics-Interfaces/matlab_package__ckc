function [CoV_ISI,ISI] = compute_CoV_ISI(MUPulses,fs)
%COMPUTE_COV_ISI Compute the Coefficient of Variation for ISIs
%
% Syntax:
%   CoV_ISI = ckc.compute_CoV_ISI(MUPulses);
%   CoV_ISI = ckc.compute_CoV_ISI(MUPulses,fs);
if nargin < 2
    fs = 2000;
end
if iscell(MUPulses)
    ISI = cell(numel(MUPulses),1);
    CoV_ISI = nan(numel(MUPulses),1);
    for ii = 1:numel(MUPulses)
        ISI{ii} = diff(MUPulses{ii}./fs);
        CoV_ISI(ii) = std(ISI{ii})/mean(ISI{ii});
    end
else
    ISI = diff(MUPulses./fs);
    CoV_ISI = std(ISI)/mean(ISI);
end
end