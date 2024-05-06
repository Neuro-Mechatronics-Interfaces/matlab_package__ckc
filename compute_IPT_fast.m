function [IPT, newPulses, residuals] = compute_IPT_fast(uni, proj, coeff, options)
%COMPUTE_IPT_FAST  Compute impulse train for CKC update (fast version)
%
% Syntax:
%   [IPT, newPulses, residuals] = ckc.compute_IPT(uni, proj, coeff, pulses, 'Name', value, ...);
%
% Inputs:
%     uni       - nChannels x nTimesteps unipolar grid data array
%     proj      - Projection matrix (nPT x nChannels)
%     coeff     - Coefficients (nChannels x nPT) of eigenvector
%     pulses    - Pulse Instants (samples)
%
% Options:
%     'SigmaThreshold' (1,1) {mustBePositive} = 3.5;
%     'IPTThreshold' (1,1) {mustBePositive} = 0.5;
%     'MaxThresholdAdjustments' (1,1) {mustBePositive, mustBeInteger} = 4;
%     'ThresholdAdjustmentScalar' (1,1) {mustBePositive, mustBeInRange(options.ThresholdAdjustmentScalar, 0.25, 0.95)} = 0.75;
%
% Output
%   IPT         - Impulse train
%   newPulses   - New predicted pulse times
%   newThreshold - Updates threshold in case insufficient suprathreshold pulses.
%
% See also: Contents

arguments
    uni
    proj
    coeff
    options.SigmaThreshold (1,1) {mustBePositive} = 3.5;
    options.IPTThreshold (1,:) {mustBePositive} = 0.5;
    options.MaxThresholdAdjustments (1,1) {mustBePositive, mustBeInteger} = 4;
    options.ThresholdAdjustmentScalar (1,1) {mustBePositive, mustBeInRange(options.ThresholdAdjustmentScalar, 0.25, 0.95)} = 0.75;
end

IPT = proj * uni;
sigma = options.SigmaThreshold.*std(IPT,[],2);
warning('off','signal:findpeaks:largeMinPeakHeight');
newThreshold = zeros(1,size(IPT,1));
newPulses = cell(1,size(IPT,1));
residuals = uni;
for ii = 1:size(IPT,1)
    IPT(ii,:) = proj(ii,:) * residuals;
    mask = IPT(ii,:) < sigma(ii);
    IPT(ii,mask) = zeros(1,sum(mask));
    IPT(ii,:) = IPT(ii,:) ./ max(abs(IPT(ii,:)));
    iCounter = 0;
    newPulses{ii} = [];
    newThreshold(ii) = options.IPTThreshold(ii);
    while (iCounter < options.MaxThresholdAdjustments)
        [~,newPulses{ii}] = findpeaks(IPT(ii,:),'MinPeakHeight',newThreshold(ii));
        if ~isempty(newPulses{ii})
            break;
        else
            newThreshold(ii) = newThreshold(ii) * options.ThresholdAdjustmentScalar;
        end
    end
    residuals = residuals - mean(residuals(:,newPulses{ii}),2).*coeff * IPT;
end
warning('on','signal:findpeaks:largeMinPeakHeight');
residuals = uni - coeff * IPT;
end