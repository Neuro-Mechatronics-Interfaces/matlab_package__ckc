function data = reader_256_poly5(filepath,filename,epoch_length,options)
%READER_256_POLY5 For reading surface EMG acquired by TMSi-SAGA using Polybench and 4x 64 8x8 standard PCB grids.
% 
% SYNTAX:
%   data = ckc.reader_256_poly5(filepath, filename, epoch_length, 'Name', value, ...);
%
% INPUTS:
%   - filepath: directory with the SIG file to be loaded.
%   - filename: EMG file to be loaded.
%   - epoch_length: (optional) length of the epoch of signal to be loaded (in s)
%
% OPTIONS:
%   ChannelMap {mustBeMember(options.ChannelMap,{'default_tmsi_pcb_64.mat','default_tmsi_pcb_128.mat','default_tmsi_pcb_256.mat'})} = 'default_tmsi_pcb_256.mat';
%   RawChannelMap = [];
%   ID = [];
%   IED (1,1) {mustBeNumeric, mustBePositive} = 8; 
%   Montage {mustBeTextScalar, mustBeMember(options.Montage,{'MONO','SD'})} = 'SD'; 
%   SampleRate (1,1) {mustBeNumeric, mustBePositive} = 4000;
%   Description {mustBeTextScalar} = "No description given.";
%   AUXchannels = [];
%   REFchannels = [];
%
% OUTPUT:
%   data: structure with the following fields;
%       SIG - two dimensional cell array with surface EMG channel in each
%             cell - SIG{r,c} is the channel in row r and column c. Missing
%             electrodes are denoted by empty arrays, e.g. SIG{1,1} = [];
%       fsamp - sampling frequency of sEMG
%       signal_length - length of a surface EMG signals (in samples)
%       montage - montage of electrodes - 'MONO' for monopolar, 'SD' for
%                 single differential
%       IED - inter-electrode distance (in mm)
%       ref_signal - measured reference signal (e.g. force) when avalable, empty array otherwise
%       AUXchannels - auxilary channels (currently not used by DEMUSEtool
%       description - text describing the data
%
% Adapted: Max Murphy 2024-04-25

arguments
    filepath
    filename
    epoch_length (1,1) = inf;
    options.ChannelMap {mustBeMember(options.ChannelMap,{'default_tmsi_pcb_64.mat','default_tmsi_pcb_128.mat','default_tmsi_pcb_256.mat'})} = 'default_tmsi_pcb_256.mat';
    options.RawChannelMap = [];
    options.ID = [];
    options.IED (1,1) {mustBeNumeric, mustBePositive} = 8; 
    options.Montage {mustBeTextScalar, mustBeMember(options.Montage,{'MONO','SD'})} = 'SD'; 
    options.SampleRate (1,1) {mustBeNumeric, mustBePositive} = 4000;
    options.Description {mustBeTextScalar} = "No description given.";
    options.AUXchannels = [];
    options.REFchannels = [];
end

f=load(fullfile(filepath,filename));

if ~isinf(epoch_length)
    f.uni = f.uni(:,1:min(size(f.uni,2),epoch_length));
end

if ~isempty(options.RawChannelMap)
    channelMap = options.RawChannelMap;
else
    channelMap = getfield(load(options.ChannelMap,'channelMap'),'channelMap');
end

data.SIG = cell(size(channelMap));
n = size(f.uni,2);
for ii=1:numel(i_elec) % electrode array
    if channelMap(ii) > 0
        data.SIG{ii} = f.uni(channelMap(ii),:);   
    else
        data.SIG{ii} = zeros(1,n);
    end
end
data.signal_length = n;monk

if isempty(options.REFchannels)
    if isfield(f, 'sync')
        if options.array_id > 2
            data.ref_signal = f.sync(2,:);
        else
            data.ref_signal = f.sync(1,:);
        end
    else
        data.ref_signal = [];
    end
else
    data.ref_signal = options.REFchannels;
end

if isfield(f,'montage')
    if ismember(f.montage,{'MONO','SD'})
        data.montage = f.montage;
    else
        data.montage = options.Montage;
    end
else
    data.montage = options.Montage;
end

if isfield(f,'grid_spacing')
    data.IED = f.grid_spacing;
else
    data.IED = options.IED;
end

if isfield(f,'sample_rate')
    data.fsamp = f.sample_rate;
else
    data.fsamp = options.SampleRate;
end

if isempty(options.AUXchannels)
    if isfield(f, 'x') && isfield(f,'y')
        data.AUXchannels = [f.x; f.y];
    else
        data.AUXchannels = [];
    end
else
    data.AUXchannels = options.AUXchannels;
end

if isfield(f,'description')
    data.description = f.description;
else
    data.description = options.Description;
end
end