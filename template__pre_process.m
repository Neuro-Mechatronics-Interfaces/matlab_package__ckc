function template__pre_process(output_filename, input_filename, options)
%TEMPLATE__PRE_PROCESS Template for pre-processing poly5 recordings, so they are ready for CKC reader.
%
% Syntax:
%   ckc.template__pre_process(output_filename, input_filename, 'Name', value, ...);
%
% Inputs:
%   output_filename (1,1) string    Output filename to save. 
%   input_filename (:,1) string     Array of strings that are input files to read/combine
%
% Options:
%   DataRoot {mustBeTextScalar} = "";
%   ApplyFilter (1,1) logical = true;
%   HighpassFilterCutoff (1,1) double = 100;
%   ApplyRMSCutoff (1,1) logical = false;
%   RMSCutoff (1,2) double = [1, 100];
%   ApplyGridInterpolation (1,1) logical = true;
%   InitialPulseOffset (1,1) {mustBeInteger} = 0; % Samples prior to first rising pulse, to include.
%   SampleRate (1,1) double {mustBeMember(options.SampleRate, [2000, 4000])} = 2000;
%   TriggerChannelIndicator {mustBeTextScalar} = 'TRIG';
%   TriggerBitMask = [];
%   ExcludedPulseIndices (1,:) {mustBeInteger,mustBePositive} = [];
%   
%
% Output:
%   Saves the concatenated and synchronized monopolar array data from all
%   input files to the file indicated by `output_filename,` in an
%   appropriate which should be compatible with the corresponding 
%   `ckc.reader_<pipeline>.m` file.  
%
% See also: Contents

arguments
    output_filename (1,1) string
    input_filename (:,1) string
    options.DataRoot {mustBeTextScalar} = "";
    options.ApplyFilter (1,1) logical = true;
    options.HighpassFilterCutoff (1,1) double = 100;
    options.ApplyRMSCutoff (1,1) logical = false;
    options.RMSCutoff (1,2) double = [1, 100];
    options.ApplyGridInterpolation (1,1) logical = true;
    options.InitialPulseOffset (1,1) {mustBeInteger} = 0; % Samples prior to first rising pulse, to include.
    options.SampleRate (1,1) double {mustBeMember(options.SampleRate, [2000, 4000])} = 2000;
    options.TriggerChannelIndicator {mustBeTextScalar} = 'TRIG';
    options.TriggerBitMask = [];
    options.ExcludedPulseIndices (1,:) {mustBeInteger,mustBePositive} = [];
    options.Description {mustBeTextScalar} = "Template function for CKC pre-processing."
end

if exist(output_filename,'file')~=0
    result = questdlg(sprintf('Output file (%s) already exists. Overwrite existing file?', output_filename), ...
        'Overwrite Existing Pre-Processing?', ...
        'Yes', 'No', 'No');
    if ~strcmpi(result,'Yes')
        disp("Pre-processing canceled.");
        return;
    end
end
if strlength(options.DataRoot) == 0
    poly5_files = strings(numel(input_filename),1);
    for ii = 1:numel(poly5_files)
        poly5_files(ii) = fullfile(options.DataRoot, input_filename(ii));
    end
else
    poly5_files = input_filename;
end

data = io.load_align_saga_data_many(poly5_files, ...
    'ApplyFilter', options.ApplyFilter, ...
    'HighpassFilterCutoff', options.HighpassFilterCutoff, ...
    'ApplyRMSCutoff', options.ApplyRMSCutoff, ...
    'RMSCutoff', options.RMSCutoff, ...
    'ApplyGridInterpolation', options.ApplyGridInterpolation, ...
    'ApplySpatialLaplacian', false, ...
    'InitialPulseOffset', options.InitialPulseOffset, ...
    'SampleRate', options.SampleRate, ...
    'TriggerChannelIndicator', options.TriggerChannelIndicator, ...
    'TriggerBitMask', options.TriggerBitMask, ...
    'ExcludedPulseIndices', options.ExcludedPulseIndices);
[iUni,iBip,iTrig] = ckc.get_saga_channel_masks(data.channels,...
    'ReturnNumericIndices',true);
uni = data.samples(iUni,:);
sample_rate = data.sample_rate;
ii = 1;
sync = data.samples(iTrig(1),:);
while ((ii < numel(iTrig)) && (numel(unique(sync))<2))
    ii = ii + 1;
    sync = data.samples(iTrig(ii),:);
end
if numel(iBip) > 0
    aux = data.samples(iBip(1),:);
else
    aux = [];
end
description = options.Description;
[p,~,~] = fileparts(output_filename);
if strlength(p) > 0
    if exist(p,'dir')==0
        mkdir(p);
    end
end
save(output_filename,'uni','aux','sync','sample_rate','description','-v7.3');

end