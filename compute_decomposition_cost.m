function [cost, PNR, CoV, L2, K, Sigma, PulseNoise] = compute_decomposition_cost(uni, MUPulses, IPTs, options)
%COMPUTE_DECOMPOSITION_COST  Cost function for CKC decomposition
%
% Syntax:
%   [cost, PNR, CoV, L2, K, Sigma, PulseNoise] = ckc.compute_decomposition_cost(uni, MUPulses, IPTs, 'Name', value, ...);
%
% Inputs:
%   uni      - Original samples matrix
%   MUpulses - Pulse sample instants
%   IPTs     - MUAP spike trains
%   
% Options:
%   'WeightCoV' - Default 1.0
%   'WeightPNR' - Default 1.0
%   'WeightKurtosis' - Default 0.1
%   'L2NormWeight' - Default 1.0
%   'SampleRate' - Default 4000
%
% Output:
%   cost - 1 x nMUAPs vector of costs
%   PNR - Pulse-to-noise ratio (dB)
%   CoV - Coefficient of variation for MUAPs ISIs
%   L2  - The L2 norm regularizer penalty
%   K        - Kurtosis (counts zeroed samples)
%   Sigma    - Variance in peaks (does not count zeroed samples)

arguments
    uni
    MUPulses
    IPTs
    options.WeightCoV (1,1) double = 1.0;
    options.WeightPNR (1,1) double = 1.0;
    options.WeightPulseNoise (1,1) double = 0.1;
    options.WeightKurtosis (1,1) double = 0.1;
    options.WeightVariance (1,1) double = 0.1;
    options.L2NormWeight (1,1) double = 1.0;
    options.NoiseUpperBound (1,1) double = 0.25;
    options.SampleRate (1,1) double = 4000;
end
n = size(IPTs,1);
cost = nan(1,n);
CoV = nan(1,n);
PNR = nan(1,n);
L2 = nan(1,n);
K = nan(1,n);
Sigma = nan(1,n);
PulseNoise = nan(1,n);
for ii = 1:n
    PNR(ii) = ckc.compute_PNR(MUPulses{ii},IPTs(ii,:),options.SampleRate);
    CoV(ii) = ckc.compute_CoV_ISI(MUPulses{ii});
    L2(ii) = norm(mean(uni(:,MUPulses{ii}),2));
    K(ii) = kurtosis(IPTs(ii,:),1);
    nonZeroMask = abs(IPTs(ii,:))>0;
    Sigma(ii) = var(IPTs(ii,nonZeroMask));
    noiseMask = IPTs(ii,:) < options.NoiseUpperBound;
    PulseNoise(ii) = sum(abs(IPTs(ii,noiseMask)));
    cost(ii) =  options.L2NormWeight*L2(ii) + options.WeightCoV * CoV(ii) - options.WeightPNR * PNR(ii) - options.WeightKurtosis*K(ii) + options.WeightVariance*Sigma(ii) + options.WeightPulseNoise * PulseNoise(ii);
end

end