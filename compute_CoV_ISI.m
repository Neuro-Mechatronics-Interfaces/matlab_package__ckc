function CoV_ISI = compute_CoV_ISI(MUPulses)
%COMPUTE_COV_ISI Compute the Coefficient of Variation for ISIs
%
% Syntax:
%   CoV_ISI = ckc.compute_CoV_ISI(MUPulses);

ISI = diff(MUPulses);
CoV_ISI = std(ISI)/mean(ISI);
end