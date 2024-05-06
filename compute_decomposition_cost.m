function [cost, PNR, CoV, L2, K, Sigma, PulseNoise, RatePenalty] = compute_decomposition_cost(uni, MUPulses, IPTs, fsamp, options)
%COMPUTE_DECOMPOSITION_COST  Cost function for CKC decomposition
%
% Syntax:
%   [cost, PNR, CoV, L2, K, Sigma, PulseNoise, RatePenalty] = ckc.compute_decomposition_cost(uni, MUPulses, IPTs, fsamp, 'Name', value, ...);
%
% Inputs:
%   uni      - Original samples matrix
%   MUpulses - Pulse sample instants
%   IPTs     - MUAP spike trains
%   fsamp    - Sample rate
%
% Options:
%   'WeightCoV' (1,1) double = 0.25;
%   'WeightPNR' (1,1) double = -0.05;
%   'WeightPulseNoise' (1,1) double = 0.01;
%   'WeightKurtosis' (1,1) double = -0.005;
%   'WeightVariance' (1,1) double = 20.0;
%   'WeightRatePenalty' (1,1) double = 0.25;
%   'WeightL2Norm' (1,1) double = 0.005;
%   'NoiseUpperBound' (1,1) double = 0.25;
%   'MaxNonPenalizedRate' (1,1) double = 50;
%
% Output:
%   cost        - 1 x nMUAPs vector of costs
%   PNR         - Pulse-to-noise ratio (dB)
%   CoV         - Coefficient of variation for MUAPs ISIs
%   L2          - The L2 norm regularizer penalty
%   K           - Kurtosis (counts zeroed samples)
%   Sigma       - Variance in peaks (does not count zeroed samples)
%   RatePenalty - Penalty for over-high rates

arguments
    uni
    MUPulses
    IPTs
    fsamp
    options.WeightCoV (1,1) double = 0.25;
    options.WeightPNR (1,1) double = -0.05;
    options.WeightPulseNoise (1,1) double = 0.01;
    options.WeightKurtosis (1,1) double = -0.005;
    options.WeightVariance (1,1) double = 20.0;
    options.WeightRatePenalty (1,1) double = 0.25;
    options.WeightL2Norm (1,1) double = 0.005;
    options.NoiseUpperBound (1,1) double = 0.25;
    options.MaxNonPenalizedRate (1,1) double = 50;
end
n = size(IPTs,1);
cost = nan(1,n);
CoV = nan(1,n);
PNR = nan(1,n);
L2 = nan(1,n);
K = nan(1,n);
Sigma = nan(1,n);
PulseNoise = nan(1,n);
RatePenalty = nan(1,n);
for ii = 1:n
    PNR(ii) = ckc.compute_PNR(MUPulses{ii},IPTs(ii,:),fsamp);
    CoV(ii) = ckc.compute_CoV_ISI(MUPulses{ii});
    L2(ii) = norm(mean(uni(:,MUPulses{ii}),2));
    K(ii) = kurtosis(IPTs(ii,:),1);
    nonZeroMask = abs(IPTs(ii,:))>0;
    Sigma(ii) = var(IPTs(ii,nonZeroMask));
    noiseMask = IPTs(ii,:) < options.NoiseUpperBound;
    PulseNoise(ii) = sum(abs(IPTs(ii,noiseMask)));
    IDR = fsamp./diff(MUPulses{ii});
    RatePenalty(ii) = sum(IDR > options.MaxNonPenalizedRate);
    cost(ii) =  options.WeightL2Norm*L2(ii) + options.WeightCoV * CoV(ii) + options.WeightPNR * PNR(ii) + options.WeightKurtosis*K(ii) + options.WeightVariance*Sigma(ii) + options.WeightPulseNoise * PulseNoise(ii) + options.WeightRatePenalty * RatePenalty(ii);
end

end