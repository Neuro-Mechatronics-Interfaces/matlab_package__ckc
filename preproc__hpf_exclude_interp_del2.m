function uni_f = preproc__hpf_exclude_interp_del2(uni, options)
%PREPROC__HPF_EXCLUDE_INTERP_DEL2 Preprocess using time-difference HPF, RMS-based exclusion, 2D interpolation, and discrete Laplacian.
%
% Example:
%   data = TMSiSAGA.Poly5.read('Max_2024_03_30_B_22.poly5');
%   uni = ckc.preproc__hpf_exclude_interp_del2(data.samples(2:65,:));   
%
% Syntax:
%   uni_f = ckc.preproc__hpf_exclude_interp_del2(uni, 'Name', value, ...);
%
% Inputs:
%   uni - Unipolar grid data 
%   
% Options:
%   'GridSize' (1,2) {mustBePositive, mustBeInteger} = [8 8];
%   'InterpolationMethod' {mustBeMember(options.InterpolationMethod, ["nearest", "linear", "v4", "natural", "cubic", "movmean", "movmedian"])} = "linear";
%
% Output:
%   uni_f - Filtered unipolar grid data
%
% See also: Contents

arguments
    uni (:, :) {mustBeNumeric}
    options.GridSize (1,2) {mustBePositive, mustBeInteger} = [8 8];
    options.InterpolationMethod {mustBeMember(options.InterpolationMethod, ["nearest", "linear", "v4", "natural", "cubic", "movmean", "movmedian"])} = "linear";
end
n_samples = size(uni, 2);
uni_f = gradient(uni); % "HPF"
uni_f(rms(uni_f,2) < 1, :) = missing;     % Indicate outlier channels
uni_f = reshape(uni_f, options.GridSize(1), options.GridSize(2), n_samples); % Put in grid form
for ii = 1:n_samples
    uni_f(:,:,ii) = fillmissing2(uni_f(:,:,ii),options.InterpolationMethod);
end
uni_f = reshape(uni_f,options.GridSize(1)*options.GridSize(2), n_samples);
end