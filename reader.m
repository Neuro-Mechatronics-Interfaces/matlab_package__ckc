function data = reader(filepath,filename,~,varargin)
%READER reads surface EMG acquired by TMSi-SAGA and arbitrary textile (monopolar 32-channel) surface array.
% 
% SYNTAX:
%   data = ckc.reader(filepath, filename, epoch_length, 'Name', value, ...);
%
% INPUTS:
%   - filepath: directory with the SIG file to be loaded.
%   - filename: EMG file to be loaded.
%   - varargin: (Optional) 'Name',value input argument pairs
%       - 'array_id': Should be 1 | 2 | 3 | 4 (numeric scalar) indicating 
%               which textile is to be read. This will only read/use data 
%               from 32 channels at a time, corresponding to one of the 
%               textile arrays from any given recording. DEFAULT: 3 
%                   (distal forelimb extrinsic wrist extensors)
%       - 'montage'     : default is 'MONO'
%       - 'IED'         : default is 8 (mm)
%       - 'fsamp'       : default is 4000 (samples/sec)
%       - AUXChannels   : default is [] (can be some ref vector; if left 
%                           empty then this becomes the [x;y] potentiometer 
%                           data tracking wrist movements)
%       - ref_signal    : default is [] (can be some ref vector; if left
%                           empty then this becomes the [sync] vector
%                           tracking task state) -- in this case a 5-bit
%                           signal tracking 
%                           (MSB)                                   (LSB)
%                               MOVE | O-HIT | O-VIS | I-HIT | I-VIS
%                       Where MOVE is move onset; O-HIT is outer-target
%                       HIT, O-VIS is outer-target visible; I-HIT is
%                       inner-target HIT; I-VIS is inner-target visible.
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
% -----------------------------------------------------------------------
% Author: Ales Holobar (ales.holobar@um.si)
% Last modified: 9. 5. 2022
%
% Adapted: Max Murphy 2022-12-13

p = inputParser;
p.addParameter('array_id', 3);
p.addParameter('montage', 'MONO');
p.addParameter('IED', 8);
p.addParameter('fsamp', 4000);
p.addParameter('AUXchannels', []);
p.addParameter('description', 'No description given.');
p.addParameter('ref_signal', []);
p.parse(varargin{:});

f=load([filepath filename]); 

% Textile array indexing:
[grid_crds(:,1), grid_crds(:,2)] = ind2sub([8,4], 1:32);
i_elec = (1:32) + ((p.Results.array_id-1) * 32);

data.SIG = {};
for ii=1:numel(i_elec) % electrode array     
  data.SIG{grid_crds(ii,1), grid_crds(ii,2)} = f.uni(i_elec(ii),:);      
end
data.signal_length = size(f.uni, 2);

if isempty(p.Results.ref_signal)
    if p.Results.array_id > 2
        data.ref_signal = f.sync(2,:);
    else
        data.ref_signal = f.sync(1,:);
    end
else
    data.ref_signal = p.Results.ref_signal;
end

data.montage = p.Results.montage;
data.IED = p.Results.IED;
data.fsamp = p.Results.fsamp;
if isempty(p.Results.AUXchannels)
    data.AUXchannels = [f.x; f.y];
else
    data.AUXchannels = p.Results.AUXchannels;
end
data.description = p.Results.description;

end