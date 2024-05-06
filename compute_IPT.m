function [IPT, newPulses, newThreshold] = compute_IPT(Ye, iRy, pulses, options)
%COMPUTE_IPT  Compute impulse train for CKC update
%
% Syntax:
%   [IPT, newPulses, newThreshold] = ckc.compute_IPT(Ye, iRy, pulses, 'Name', value, ...);
%
% Inputs:
%     Ye        - Extended signal vector
%     iRy       - (Pseudo-)Inverse covariance matrix
%     pulses    - Pulse Instants (samples)
%
% Options:
%     'ExtensionFactor' (1,1) {mustBePositive, mustBeInteger} = 40;
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
    Ye
    iRy
    pulses
    options.ExtensionFactor (1,1) {mustBePositive, mustBeInteger} = 40;
    options.SigmaThreshold (1,1) {mustBePositive} = 3.5;
    options.IPTThreshold (1,1) {mustBePositive} = 0.5;
    options.MaxThresholdAdjustments (1,1) {mustBePositive, mustBeInteger} = 4;
    options.ThresholdAdjustmentScalar (1,1) {mustBePositive, mustBeInRange(options.ThresholdAdjustmentScalar, 0.25, 0.95)} = 0.75;
end

IPT = sum(Ye(:,pulses),2)'*iRy*Ye;
IPT([1:(2*options.ExtensionFactor),(end-(2*options.ExtensionFactor)):end]) = 0;
sigma = options.SigmaThreshold.*std(IPT);
mask = IPT < sigma;
IPT(mask) = zeros(1,sum(mask));
IPT = IPT ./ max(abs(IPT));
warning('off','signal:findpeaks:largeMinPeakHeight');
iCounter = 0;
newPulses = [];
newThreshold = options.IPTThreshold;
while (iCounter < options.MaxThresholdAdjustments)
    [~,newPulses] = findpeaks(IPT,'MinPeakHeight',newThreshold);
    if ~isempty(newPulses)
        break;
    else
        newThreshold = newThreshold * options.ThresholdAdjustmentScalar;
    end
end
warning('on','signal:findpeaks:largeMinPeakHeight');
end