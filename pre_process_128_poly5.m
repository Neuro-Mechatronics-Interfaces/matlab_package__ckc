function pre_process_128_poly5(output_filename, poly5_file_extensor, poly5_file_flexor, options)
%PRE_PROCESS_128_POLY5  For pre-processing sets of 4x poly5 recordings, using 8x8 standard PCB grids and acquired using polybench, so that they are ready for the CKC reader.
%
% Syntax:
%   ckc.pre_process_128_poly5(output_filename, poly5_file_extensor, poly5_file_flexor, 'Name', value, ...);
%
% Inputs:
%   output_filename (1,1) string
%   poly5_file_extensor (1,1) string
%   poly5_file_flexor (1,1) string
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
    poly5_file_extensor (1,1) string
    poly5_file_flexor (1,1) string
    options.DataRoot {mustBeTextScalar} = "";
    options.AlignSync (1,1) logical = true;
    options.ApplyFilter (1,1) logical = true;
    options.ApplyCAR (1,1) logical = true;
    options.ApplySpatialLaplacian (1,1) logical = false;
    options.HighpassFilterCutoff (1,1) double = 100;
    options.ApplyRMSCutoff (1,1) logical = true;
    options.RMSCutoff (1,2) double = [1, 100];
    options.ZeroMissing (1,1) logical = true; % Sets "missing" samples as zeros
    options.ApplyGridInterpolation (1,1) logical = true;
    options.InitialPulseOffset (1,1) {mustBeInteger} = 0; % Samples prior to first rising pulse, to include.
    options.SampleRate (1,1) double {mustBeMember(options.SampleRate, [2000, 4000])} = 2000;
    options.Sync = []; % If specified, should be a 2 x k array where first row is time and second is sync value
    options.SyncTarget = [];
    options.InvertSyncLogic = [];
    options.TriggerChannelIndicator {mustBeTextScalar} = 'TRIG';
    options.TriggerBitMask = 1;
    options.ExcludedPulseIndices (1,:) {mustBeInteger,mustBePositive} = [];
    options.IsTextile64 (1,1) logical = true;
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
    poly5_files = [poly5_file_extensor; ...
                   poly5_file_flexor];
else
    poly5_files = [fullfile(options.DataRoot,poly5_file_extensor); ...
                   fullfile(options.DataRoot,poly5_file_flexor)];
end
if isscalar(options.TriggerBitMask)
    triggerBitMask = ones(numel(poly5_files),1).*options.TriggerBitMask;
else
    triggerBitMask = options.TriggerBitMask;
end
[data,~,ch_name] = io.load_align_saga_data_many(poly5_files, ...
    'ApplyFilter', options.ApplyFilter, ...
    'ApplyCAR', options.ApplyCAR, ...
    'HighpassFilterCutoff', options.HighpassFilterCutoff, ...
    'ApplyRMSCutoff', options.ApplyRMSCutoff, ...
    'RMSCutoff', options.RMSCutoff, ...
    'ZeroMissing',options.ZeroMissing,...
    'ApplyGridInterpolation', options.ApplyGridInterpolation, ...
    'ApplySpatialLaplacian', options.ApplySpatialLaplacian, ...
    'InitialPulseOffset', options.InitialPulseOffset, ...
    'InvertLogic', options.InvertSyncLogic, ...
    'SampleRate', options.SampleRate, ...
    'TriggerChannelIndicator', options.TriggerChannelIndicator, ...
    'TriggerBitMask', triggerBitMask, ...
    'IsTextile64', options.IsTextile64, ...
    'ExcludedPulseIndices', options.ExcludedPulseIndices);
[iUni,~,iTrig] = ckc.get_saga_channel_masks(ch_name,...
    'ReturnNumericIndices',true);
uni = data.samples(iUni,:);
sample_rate = data.sample_rate;
ii = 1;
sync = data.samples(iTrig(1),:);
if isempty(options.InvertSyncLogic)
    sync = double(bitand(data.samples,triggerBitMask(1))==triggerBitMask(1));
else
    if numel(options.InvertSyncLogic) > 1

    else
        if options.InvertSyncLogic
            sync = double(bitand(data.samples,triggerBitMask(1))~=triggerBitMask(1));
        else
            sync = double(bitand(data.samples,triggerBitMask(1))==triggerBitMask(1));
        end
    end
end

while ((ii < numel(iTrig)) && (numel(unique(sync))<2))
    ii = ii + 1;
    if isempty(options.InvertSyncLogic)
        sync = double(bitand(data.samples,triggerBitMask(ii))==triggerBitMask(ii));
    else
        if numel(options.InvertSyncLogic) > 1
    
        else
            if options.InvertSyncLogic
                sync = double(bitand(data.samples,triggerBitMask(ii))~=triggerBitMask(ii));
            else
                sync = double(bitand(data.samples,triggerBitMask(ii))==triggerBitMask(ii));
            end
        end
    end
end
t_data = 0:(1/data.sample_rate):((numel(sync)-1)/data.sample_rate);
t_start = t_data(find(sync > 0,1,'first'));
if ~isempty(options.Sync)
    if options.AlignSync
        sync_in = [0, 0, options.Sync(2,:), 0, 0];
        sync_time = options.Sync(1,:) + t_start;
        sync_time = [0, sync_time(1)-(1/data.sample_rate), sync_time, sync_time(end)+(1/data.sample_rate), t_data(end)];
        sync = interp1(sync_time, sync_in, t_data, 'linear');
    else
        sync = options.Sync;
    end
end
if ~isempty(options.SyncTarget)
    sync_in = [0, options.SyncTarget(:,2)', 0];
    sync_time = options.SyncTarget(:,1)' + t_start;
    [sync_time, ikeep] = unique([0, sync_time, t_data(end)]);
    sync_out = interp1(sync_time, sync_in(ikeep), t_data, 'previous');
    sync = [sync; sync_out];
end

aux = [];
description = options.Description;
[p,~,~] = fileparts(output_filename);
if strlength(p) > 0
    if exist(p,'dir')==0
        mkdir(p);
    end
end
save(output_filename,'uni','aux','sync','sample_rate','description','-v7.3');

end