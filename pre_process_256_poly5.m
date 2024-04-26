function pre_process_256_poly5(output_filename, poly5_file_proximal_extensor, poly5_file_distal_extensor, poly5_file_proximal_flexor, poly5_file_distal_flexor, options)
%PRE_PROCESS_256_POLY5  For pre-processing sets of 4x poly5 recordings, using 8x8 standard PCB grids and acquired using polybench, so that they are ready for the CKC reader.
%
% Syntax:
%   ckc.template__pre_process(output_filename, poly5_file_proximal_extensor, poly5_file_distal_extensor, poly5_file_proximal_flexor, poly5_file_distal_flexor, 'Name', value, ...);
%
% Inputs:
%   output_filename (1,1) string
%   poly5_file_proximal_extensor (1,1) string
%   poly5_file_distal_extensor (1,1) string
%   poly5_file_proximal_flexor (1,1) string
%   poly5_file_distal_flexor (1,1) string
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
    poly5_file_proximal_extensor (1,1) string
    poly5_file_distal_extensor (1,1) string
    poly5_file_proximal_flexor (1,1) string
    poly5_file_distal_flexor (1,1) string
    options.DataRoot {mustBeTextScalar} = "";
    options.ApplyFilter (1,1) logical = true;
    options.HighpassFilterCutoff (1,1) double = 100;
    options.ApplyRMSCutoff (1,1) logical = false;
    options.RMSCutoff (1,2) double = [1, 100];
    options.ApplyGridInterpolation (1,1) logical = true;
    options.InitialPulseOffset (1,1) {mustBeInteger} = 0; % Samples prior to first rising pulse, to include.
    options.SampleRate (1,1) double {mustBeMember(options.SampleRate, [2000, 4000])} = 2000;
    options.TriggerChannelIndicator {mustBeTextScalar} = 'TRIG';
    options.TriggerBitMask = 1;
    options.ExcludedPulseIndices (1,:) {mustBeInteger,mustBePositive} = [];
    options.Description {mustBeTextScalar} = "1-64 = PROX-EXT; 65-128 = DIST-EXT; 129-192 = PROX-FLX; 193-256 = DIST-FLX."
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
    poly5_files = [poly5_file_proximal_extensor; poly5_file_distal_extensor; poly5_file_proximal_flexor; poly5_file_distal_flexor];
else
    poly5_files = [fullfile(options.DataRoot,poly5_file_proximal_extensor); ...
                   fullfile(options.DataRoot,poly5_file_distal_extensor); ...
                   fullfile(options.DataRoot,poly5_file_proximal_flexor); ...
                   fullfile(options.DataRoot,poly5_file_distal_flexor)];
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
