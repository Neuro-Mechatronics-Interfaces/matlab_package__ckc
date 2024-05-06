function [coeff, proj, explained] = uni_2_pcs(uni, options)
%UNI_2_PCS  Return PCs from unipolar array
%
% Syntax:
%   [coeff, proj, explained] = ckc.uni_2_pcs(uni, "Name", value, ...);
%
% Inputs:
%   uni - 64 x nSamples data array
%  
% Options:
%   MinPCVarCapt (1,1) {mustBePositive, mustBeInRange(0, 100)} = 2
%   MinPCsKept (1,1) {mustBePositive, mustBeInteger} = 3;
%
% Output:
%   coeff - 64 x nPCs coefficients, each col takes us from IPT to uni array
%   proj  - nPCs x 64 coefficients, each row takes us from uni to
%               eigenvector projection (--> IPT)
%   explained - Total variance percent explained from these PCs

arguments
    uni (64,:)
    options.MinPCVarCapt (1,1) double {mustBePositive, mustBeInRange(options.MinPCVarCapt, 0, 100)} = 2; % Percent of data explained requirement
    options.MinPCsKept (1,1) double {mustBePositive, mustBeInteger} = 3;
end

warning('off','stats:pca:ColRankDefX');
[~, coeff, ~, ~, explained] = pca(uni);
n = find(explained > options.MinPCVarCapt, 1, 'last');
if isempty(n)
    n = options.MinPCsKept;
else
    n = max(n, options.MinPCsKept);
end
coeff = coeff(:,1:n);
explained = sum(explained(1:n));
proj = pinv(coeff);
warning('on','stats:pca:ColRankDefX');

end