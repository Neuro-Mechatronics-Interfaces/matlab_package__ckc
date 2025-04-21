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
    options.Debug (1,1) logical = false;
    options.UseFirstSampleIfNoSyncPulse (1,1) logical = true;
    options.AlignSync (1,1) logical = true;
    options.ApplyFilter (1,1) logical = true;
    options.ApplyCAR (1,1) logical = true;
    options.ApplySpatialFilter (1,1) logical = false;
    options.ChannelMap (1,256) double {mustBePositive, mustBeInteger, mustBeInRange(options.ChannelMap,1,256)} = 1:256;
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
    options.IsTextile64 (:,1) logical = true;
    options.SwappedTextileCables (1,4) logical = false(1,4);
    options.UsePulseSync (1,1) logical = false;
    options.MaxSyncPulseWidth (1,1) double = 0.6;
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
    poly5_files = [string(poly5_file_proximal_extensor); ...
                   string(poly5_file_distal_extensor); ...
                   string(poly5_file_proximal_flexor); ...
                   string(poly5_file_distal_flexor)];
else
    poly5_files = [string(sprintf('%s/%s',options.DataRoot,poly5_file_proximal_extensor)); ...
                   string(sprintf('%s/%s',options.DataRoot,poly5_file_distal_extensor)); ...
                   string(sprintf('%s/%s',options.DataRoot,poly5_file_proximal_flexor)); ...
                   string(sprintf('%s/%s',options.DataRoot,poly5_file_distal_flexor))];
end
if isscalar(options.TriggerBitMask)
    triggerBitMask = ones(numel(poly5_files),1).*options.TriggerBitMask;
else
    triggerBitMask = options.TriggerBitMask;
end
[data,~,ch_name] = io.load_align_saga_data_many(poly5_files, ...
    'ApplyFilter', options.ApplyFilter, ...
    'ApplyCAR', options.ApplyCAR, ...
    'Debug', options.Debug, ...
    'HighpassFilterCutoff', options.HighpassFilterCutoff, ...
    'ApplyRMSCutoff', options.ApplyRMSCutoff, ...
    'RMSCutoff', options.RMSCutoff, ...
    'ZeroMissing',options.ZeroMissing,...
    'ApplyGridInterpolation', options.ApplyGridInterpolation, ...
    'ApplySpatialFilter', options.ApplySpatialFilter, ...
    'InitialPulseOffset', options.InitialPulseOffset, ...
    'InvertLogic', options.InvertSyncLogic, ...
    'SampleRate', options.SampleRate, ...
    'TriggerChannelIndicator', options.TriggerChannelIndicator, ...
    'TriggerBitMask', triggerBitMask, ...
    'IsTextile64', options.IsTextile64, ...
    'SwappedTextileCables', options.SwappedTextileCables, ...
    'UseFirstSampleIfNoSyncPulse', options.UseFirstSampleIfNoSyncPulse, ...
    'UsePulseSync', options.UsePulseSync, ...
    'MaxSyncPulseWidth', options.MaxSyncPulseWidth, ...
    'ExcludedPulseIndices', options.ExcludedPulseIndices);
[iUni,iBip,iTrig] = ckc.get_saga_channel_masks(ch_name,...
    'ReturnNumericIndices',true);
uni = data.samples(iUni,:);
uni = uni(options.ChannelMap,:);

sample_rate = data.sample_rate;
ii = 1;
sync = data.samples(iTrig(1),:);
all_sync = data.samples(iTrig,:);
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
save(output_filename,'uni','aux','sync','all_sync','sample_rate','description','-v7.3');

end
