function pre_process_textiles_long_fmt(T, n_textiles, varargin)
%PRE_PROCESS_TEXTILES_LONG_FMT  Batch concatenate the textile arrays into corresponding files 
%
% Syntax:
%   ckc.pre_process_textiles_LONG_FMT(T, n_textiles, 'Name', value, ...);
%
% Inputs:
%   T - Table returned by io.read_events. It has columns flagging which
%           recording blocks are to be excluded.
%
% Output:
%   In generated_data tank folder, there will be new 
%
% Example:
%   T = io.read_events("Spencer", 2022, 12, 8);
%   ckc.pre_process_textiles_long_fmt(T, 4);
%
% See also: Contents

[generated_data_folder, raw_data_folder] = parameters('generated_data_folder', 'raw_data_folder');
p = inputParser;
p.addParameter('montage', 'MONO');
p.addParameter('IED', 8.75);
p.addParameter('fsamp', 4000);
p.addParameter('generated_data_folder', generated_data_folder);
p.addParameter('raw_data_folder', raw_data_folder);
p.addParameter('name_expr_cal', '%s_CKC-Cal_Array-%d.mat');
p.addParameter('name_expr_full', '%s_CKC-Full_Array-%d.mat');
p.parse(varargin{:});
generated_data_folder = p.Results.generated_data_folder;
raw_data_folder = p.Results.raw_data_folder;
montage = p.Results.montage;
IED = p.Results.IED;
fsamp = p.Results.fsamp;

exc = T.excluded_by_pots | T.excluded_by_noise | T.excluded_by_manual;
T = T(~exc,:);

% First, export a set of "calibration" files:
%   1/array => 2-4 total
%   1/orientation, 1/direction, 1/target => 32 total
%   => Up to 128 "calibrations"
G = findgroups(T(:, {'direction', 'target_index', 'orientation'}));
uG = unique(G);
iCal = nan(numel(uG)*4,1); % Get 4 and concatenate together for calibration files.
for k = 1:numel(uG)
    vec = (1:4) + ((k-1) * 4);
    iCal(vec) = find(G == uG(k), 4, 'first');
end

SUBJ = T.subject(1);
YYYY = year(T.date(1));
MM = month(T.date(1));
DD = day(T.date(1));
tank = sprintf('%s_%04d_%02d_%02d', SUBJ, YYYY, MM, DD);

about = io.parse_about(SUBJ, YYYY, MM, DD, 'raw_data_folder', raw_data_folder);
desc_locs = horzcat(about.SAGA.A.Location, about.SAGA.B.Location);

grid_crds = nan(32, 1);
[grid_crds(:,1), grid_crds(:,2)] = ind2sub([8,4], 1:32);

C = T(iCal,:);
nTotal = size(C,1);
fprintf(1,'Saving calibration files...000%%\n');
for ii = 1:n_textiles
    SIG = cell(8,4);
    uni = [];
    ref_signal = [];
    AUXchannels = [];
    description = strrep(strrep(char(desc_locs(ii)), ' (', ' - '), ')', '');
    i_elec = (1:32) + ((ii-1) * 32);
    i_prog = (ii/n_textiles)*100;
    for iC = 1:nTotal
        data = io.load_wrist_event_table_trial(C(iC,:), ...
            'raw_data_folder', raw_data_folder, ...
            'generated_data_folder', generated_data_folder);
        uni = horzcat(uni, data.uni(i_elec,:)); %#ok<AGROW> 
        ref_signal = horzcat(ref_signal, sqrt(double(data.x).^2 + double(data.y).^2)); %#ok<AGROW> 
        AUXchannels = horzcat(AUXchannels, [double(data.x); double(data.y); double(data.sync)]); %#ok<AGROW> 
        fprintf(1,'\b\b\b\b\b%03d%%\n', round(i_prog * (iC / nTotal)));
    end
    uni = double(uni);
    signal_length = size(uni,2);
    r = rms(uni, 2);
    iBad = find(r >= (mean(r) + 2*std(r)));
    iGood = setdiff(1:32, iBad);
    uni(iBad,:) = zeros(numel(iBad), signal_length);
    uni(iGood,:) = uni(iGood,:) - mean(uni(iGood,:), 1);
    for ik = 1:32
        SIG{grid_crds(ik,1), grid_crds(ik,2)} = uni(ik, :);
    end
    save(fullfile(generated_data_folder, SUBJ, tank,  sprintf(p.Results.name_expr_cal, tank, ii)), 'SIG', 'signal_length', 'ref_signal', 'AUXchannels', 'fsamp', 'montage', 'IED', 'description', '-v7.3');
    fprintf(1,'\b\b\b\b\b%03d%%\n', round(i_prog));
end


% Next, concatenate ALL
fprintf(1,'Saving FULL files...000%%\n');
nTotal = size(T,1);
for ii = 1:n_textiles
    SIG = cell(8,4);
    uni = [];
    ref_signal = [];
    AUXchannels = [];
    description = strrep(strrep(char(desc_locs(ii)), ' (', ' - '), ')', '');
    i_elec = (1:32) + ((ii-1) * 32);
    i_prog = (ii/n_textiles)*100;
    for iC = 1:nTotal
        data = io.load_wrist_event_table_trial(T(iC,:), ...
            'raw_data_folder', raw_data_folder, ...
            'generated_data_folder', generated_data_folder);
        uni = horzcat(uni, data.uni(i_elec,:)); %#ok<AGROW> 
        ref_signal = horzcat(ref_signal, sqrt(double(data.x).^2 + double(data.y).^2)); %#ok<AGROW> 
        AUXchannels = horzcat(AUXchannels, [double(data.x); double(data.y); double(data.sync)]); %#ok<AGROW> 
        fprintf(1,'\b\b\b\b\b%03d%%\n', round(i_prog * (iC / nTotal)));
    end
    uni = double(uni);
    signal_length = size(uni,2);
    r = rms(uni, 2);
    iBad = find(r >= (mean(r) + 2*std(r)));
    iGood = setdiff(1:32, iBad);
    uni(iBad,:) = zeros(numel(iBad), signal_length);
    uni(iGood,:) = uni(iGood,:) - mean(uni(iGood,:), 1);
    for ik = 1:32
        SIG{grid_crds(ik,1), grid_crds(ik,2)} = uni(ik, :);
    end
    save(fullfile(generated_data_folder, SUBJ, tank, sprintf(p.Results.name_expr_full, tank, ii)), 'SIG', 'signal_length', 'ref_signal', 'AUXchannels', 'fsamp', 'montage', 'IED', 'description', '-v7.3');
    fprintf(1,'\b\b\b\b\b%03d%%\n', round(i_prog));
end

end