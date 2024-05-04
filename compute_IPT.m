function [IPT, newPulses] = compute_IPT(Ye, iRy, pulses, options)
%COMPUTE_IPT  Compute impulse train for CKC update
%
% Syntax:
%   [IPT, newPulses] = ckc.compute_IPT(Ye, iRy, pulses, 'Name', value, ...);
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
%
% Output
%   IPT         - Impulse train
%   newPulses   - New predicted pulse times
%
% See also: Contents

arguments
    Ye
    iRy
    pulses
    options.ExtensionFactor (1,1) {mustBePositive, mustBeInteger} = 40;
    options.SigmaThreshold (1,1) {mustBePositive} = 3.5;
    options.IPTThreshold (1,1) {mustBePositive} = 0.5;
end

IPT = sum(Ye(:,pulses),2)'*iRy*Ye;
IPT([1:(2*options.ExtensionFactor),(end-(2*options.ExtensionFactor)):end]) = 0;
sigma = options.SigmaThreshold.*std(IPT);
mask = IPT < sigma;
IPT(mask) = zeros(1,sum(mask));
IPT = IPT ./ max(abs(IPT));
[~,newPulses] = findpeaks(IPT,'MinPeakHeight',options.IPTThreshold);

end