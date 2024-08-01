function pre_process_2x_4textile(input_root, block, triggerBitMask, options)
%PRE_PROCESS_2x_4textile  For pre-processing sets of 4x poly5 recordings, using 8x8 standard PCB grids and acquired using polybench, so that they are ready for the CKC reader.
%
% Syntax:
%   ckc.pre_process_2x_4textile(input_root, block, 'Name', value, ...);
%
% Inputs:
%   input_root (1,1) string
%   block (1,1) double {mustBeInteger}
%   triggerBitMask (1,:) double {mustBeInteger}
%
% Options:
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
    input_root {mustBeTextScalar, mustBeFolder}
    block (1,1) double {mustBeInteger}
    triggerBitMask (1,:) double {mustBeInteger} = 2;
    options.AlignSync (1,1) logical = true;
    options.ApplyFilter (1,1) logical = true;
    options.ApplyCAR (1,1) logical = true;
    options.ApplySpatialFilter (1,1) logical = false;
    options.HighpassFilterCutoff (1,1) double = 100;
    options.ApplyRMSCutoff (1,1) logical = true;
    options.RMSCutoff (1,2) double = [1, 100];
    options.ZeroMissing (1,1) logical = true; % Sets "missing" samples as zeros
    options.ApplyGridInterpolation (1,1) logical = true;
    options.InitialPulseOffset (1,1) {mustBeInteger} = 0; % Samples prior to first rising pulse, to include.
    options.SampleRate (1,1) double {mustBeMember(options.SampleRate, [2000, 4000])} = 2000;
    options.Sync = []; % If specified, should be a 2 x k array where first row is time and second is sync value
    options.SyncTarget = [];
    options.OutputFileTag = 'Synchronized';
    options.Description = "";
    options.Subfolder = "TMSi";
end

[~,TANK] = fileparts(input_root);
A = dir(fullfile(input_root, options.Subfolder, sprintf("%s_A*_%d.poly5", TANK, block)));
B = dir(fullfile(input_root, options.Subfolder, sprintf("%s_B*_%d.poly5", TANK, block)));
poly5_files = [string(fullfile(A(1).folder, A(1).name)); string(fullfile(B(1).folder, B(1).name))];
if strlength(options.Description) > 0
    description = sprintf("%s_2x_4mm_%d: %s", TANK, block, options.Description);
else
    description = sprintf("%s_2x_4mm_%d: WITH TABLET", TANK, block);
end
[data,~,ch_name] = io.load_align_saga_data_many(poly5_files, ...
    'ApplyFilter', options.ApplyFilter, ...
    'ApplyCAR', options.ApplyCAR, ...
    'HighpassFilterCutoff', options.HighpassFilterCutoff, ...
    'ApplyRMSCutoff', options.ApplyRMSCutoff, ...
    'RMSCutoff', options.RMSCutoff, ...
    'ZeroMissing',options.ZeroMissing,...
    'ApplyGridInterpolation', options.ApplyGridInterpolation, ...
    'ApplySpatialFilter', options.ApplySpatialFilter, ...
    'SampleRate', options.SampleRate, ...
    'TriggerBitMask', triggerBitMask, ...
    'IsTextile64', true);
[iUni,~,iTrig] = ckc.get_saga_channel_masks(ch_name,...
    'ReturnNumericIndices',true);
uni = data.samples(iUni,:);
sample_rate = data.sample_rate;
sync = double(bitand(data.samples(iTrig(1),:),triggerBitMask(1))==0);

aux = [];
output_filename = fullfile(input_root, options.Subfolder, sprintf("%s_%s_%d.mat", TANK, options.OutputFileTag, block));
save(output_filename,'uni','aux','sync','sample_rate','description','-v7.3');

end