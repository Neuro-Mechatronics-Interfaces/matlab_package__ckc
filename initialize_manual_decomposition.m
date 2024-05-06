function fig = initialize_manual_decomposition(uni, fsamp, options)
%INITIALIZE_MANUAL_DECOMPOSITION  Manually initializes the decomposition for CKC
%
% Example:
%   data = TMSiSAGA.Poly5.read('Max_2024_03_30_B_22.poly5');
%   uni = ckc.preproc__hpf_exclude_interp_del2(data.samples(2:65,:));   
%   [IPTs, MUPulses, info, t, R_inv] = ckc.initialize_manual_decomposition(uni, data.sample_rate);
%
% Syntax:
%   [IPTs, MUPulses, info, t, R_inv] = ckc.initialize_manual_decomposition(uni, fsamp, 'Name', value, ...);
%
% Inputs:
%   uni - Signal from grid montage EMG channels to use with CKC decomp.
%   fsamp - Sample rate of this recording
%
% Options:
%   ExtensionFactor = 40;
%   NumPCs (1,1) {mustBePositive, mustBeInteger} = 6;
%   SigmaThreshold (1,1) {mustBePositive} = 3.5;
%   IPTThreshold (1,1) {mustBePositive} = 0.5;
%   MaxThresholdAdjustments (1,1) {mustBePositive, mustBeInteger} = 4;
%   ThresholdAdjustmentScalar (1,1) {mustBePositive, mustBeInRange(options.ThresholdAdjustmentScalar, 0.25, 0.95)} = 0.75;
%   IPTThresholdStepSize (1,1) {mustBePositive} = 0.025;
%   IPTCleaningIterationsMax (1,1) {mustBePositive} = 50;
%   SubSampleVector (1,2) = [nan nan];
%   WeightCoV (1,1) double = 0.25;
%   WeightPNR (1,1) double = -0.05;
%   WeightRatePenalty (1,1) double = 0.025;
%   WeightL2Norm (1,1) double = 0.0050;;
%   WeightKurtosis (1,1) double = -0.0050;
%   WeightVariance (1,1) double = 50;
%   MaxNonPenalizedRate (1,1) double = 50;
%   MaxStrikes (1,1) {mustBeInteger, mustBePositive} = 10;
%   SampleRate (1,1) double = 4000;
%   Verbose (1,1) logical = true;

arguments
    uni (64,:) double 
    fsamp (1,1) {mustBePositive}
    % options.ExtensionFactor = 40;
    options.NumPCsMax (1,1) {mustBePositive, mustBeInteger} = 20;
    options.SigmaThreshold (1,1) {mustBePositive} = 2.5;
    options.IPTThreshold (1,1) {mustBePositive} = 0.4;
    options.MaxThresholdAdjustments (1,1) {mustBePositive, mustBeInteger} = 4;
    options.ThresholdAdjustmentScalar (1,1) {mustBePositive, mustBeInRange(options.ThresholdAdjustmentScalar, 0.25, 0.95)} = 0.75;
    options.IPTThresholdStepSize (1,1) {mustBePositive} = 0.025;
    options.IPTCleaningIterationsMax (1,1) {mustBePositive} = 50;
    options.IPTNoiseUpperBound (1,1) double = 0.25;
    options.SubSampleVector (1,2) = [nan nan];
    options.WeightCoV (1,1) double = 0.25;
    options.WeightL2Norm (1,1) double = 0.002;
    options.WeightPNR (1,1) double = -0.05;
    options.WeightPulseNoise (1,1) double = 0.01;
    options.WeightKurtosis (1,1) double = -0.005;
    options.WeightVariance (1,1) double = 50;
    options.WeightRatePenalty (1,1) double = 0.25;
    options.MaxNonPenalizedRate (1,1) double = 40;
    options.PreferredIntervalTolerance (1,1) double = 0.2;
    options.MinPCVarCapt (1,1) double {mustBePositive, mustBeInRange(options.MinPCVarCapt, 0, 100)} = 2; % Percent of data explained requirement
    options.MinPCsKept (1,1) double {mustBePositive, mustBeInteger} = 3;
    options.Verbose (1,1) logical = true;
end

fig = figure('Color', 'w', 'Name', 'Manual CKC Decomposer');
fig.UserData = struct;

fig.UserData.fsamp = fsamp;
fig.UserData.iPT = ones(2,1);
fig.UserData.iRep = 1;
fig.UserData.nRep = 1;
fig.UserData.nsamp = size(uni,2);

[fig.UserData.coeff, fig.UserData.proj] = ckc.uni_2_pcs(uni, ...
    'MinPCsKept', options.MinPCsKept, ...
    'MinPCVarCapt', options.MinPCVarCapt);
fig.UserData.nPT = ones(2,1).*size(fig.UserData.coeff,2);
fig.UserData.coeff = repmat({fig.UserData.coeff},2,1);
fig.UserData.proj = repmat({fig.UserData.proj},2,1);

fig.UserData.uni = {uni; uni};
fig.UserData.IPTs = repmat({zeros(fig.UserData.nPT(1),fig.UserData.nsamp)},2,1);
fig.UserData.IPTThreshold = repmat({ones(1,fig.UserData.nPT(1)).*options.IPTThreshold},2,1);
warning('off','signal:findpeaks:largeMinPeakHeight');


L = tiledlayout(fig, 3,6);


fig.UserData.t = (0:(fig.UserData.nsamp-1))/fig.UserData.fsamp;
fig.UserData.IPTAxes = nexttile(L,1,[2 6]);
set(fig.UserData.IPTAxes,'NextPlot','add','FontName','Tahoma','XLim',[fig.UserData.t(1), fig.UserData.t(end)],'YLim',[0 1]);
fig.UserData.IPTLine = line(fig.UserData.IPTAxes, fig.UserData.t, zeros(size(fig.UserData.t)), 'Color', 'k');
fig.UserData.MUPulseLine = line(fig.UserData.IPTAxes, nan, nan, 'LineStyle', 'none', 'Marker', 'o', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'none');
fig.UserData.IPTThresholdLine = yline(fig.UserData.IPTAxes, options.IPTThreshold, 'r:', 'IPT Threshold', ...
    'ButtonDownFcn',@ipt_threshold_line_grab_callback,'LineWidth', 2);

fig.UserData.Title = title(fig.UserData.IPTAxes, "IPT-1 | Rep-1 | N = 0 | PNR = 0 dB", 'FontName', 'Consolas', 'Color', 'k');
popupValue = strings(fig.UserData.nPT(1),1);
for k = 1:fig.UserData.nPT(1)
    popupValue(k) = sprintf("IPT-%02d",k);
end
fig.UserData.IPTPopupMenu = uicontrol(fig, ...
    'Style','popupmenu','FontName','Consolas', ...
    'Units','Normalized',...
    'Position',[0.05 0.05 0.2 0.075], ...
    'String',popupValue, ...
    'Callback',@handle_changing_ipt_index);
fig.UserData.AddRepButton = uicontrol(fig, ...
    'Style', 'pushbutton', 'FontName', 'Lock IPT', ...
    'Units','Normalized', ...
    'Position',[0.25 0.05 0.2 0.05], ...
    'String', "Add Rep", ...
    'Callback', @handle_adding_new_rep);

[~, fig.UserData.MUPulses] = ckc.uni_2_pks(...
    fig.UserData.uni{fig.UserData.iPT(fig.UserData.iRep)}, ...
    'NPCs', fig.UserData.nPT(1), ...
    'Verbose', options.Verbose);
fig.UserData.MUPulses = repmat({fig.UserData.MUPulses},2,1);
fig.UserData.Cost = repmat({zeros(1,fig.UserData.nPT(1))},2,1);
fig.UserData.PNR = repmat({zeros(1,fig.UserData.nPT(1))},2,1);
fig.UserData.CoV = repmat({zeros(1,fig.UserData.nPT(1))},2,1);
fig.UserData.L2 = repmat({zeros(1,fig.UserData.nPT(1))},2,1);
fig.UserData.Kurtosis = repmat({zeros(1,fig.UserData.nPT(1))},2,1);
fig.UserData.Variance = repmat({zeros(1,fig.UserData.nPT(1))},2,1);
fig.UserData.PulseNoise = repmat({zeros(1,fig.UserData.nPT(1))},2,1);
fig.UserData.RatePenalty = repmat({zeros(1,fig.UserData.nPT(1))},2,1);

recalculate_ipts();
 

    function recalculate_ipts(~,~)
        fig.Pointer = "watch";
        drawnow();
        [fig.UserData.IPTs{fig.UserData.iRep}, fig.UserData.MUPulses{fig.UserData.iRep},fig.UserData.uni{fig.UserData.iRep+1}] = ...
                ckc.compute_IPT_fast(fig.UserData.uni{fig.UserData.iRep}, fig.UserData.proj{fig.UserData.iRep}, fig.UserData.coeff{fig.UserData.iRep}, ...
                    'SigmaThreshold', options.SigmaThreshold, ...
                    'IPTThreshold', fig.UserData.IPTThreshold{fig.UserData.iRep}, ...
                    'MaxThresholdAdjustments', options.MaxThresholdAdjustments, ...
                    'ThresholdAdjustmentScalar', options.ThresholdAdjustmentScalar);  
        [fig.UserData.Cost{fig.UserData.iRep}, fig.UserData.PNR{fig.UserData.iRep}, fig.UserData.CoV{fig.UserData.iRep}, fig.UserData.L2{fig.UserData.iRep}, fig.UserData.Kurtosis{fig.UserData.iRep}, fig.UserData.Variance{fig.UserData.iRep}, fig.UserData.PulseNoise{fig.UserData.iRep}, fig.UserData.RatePenalty{fig.UserData.iRep}] = ...
            ckc.compute_decomposition_cost(...
                fig.UserData.uni{fig.UserData.iRep}, fig.UserData.MUPulses{fig.UserData.iRep}, fig.UserData.IPTs{fig.UserData.iRep}, fig.UserData.fsamp, ...
                'NoiseUpperBound', options.IPTNoiseUpperBound, ...
                'WeightL2Norm', options.WeightL2Norm, ...
                'WeightCoV', options.WeightCoV, ...
                'WeightKurtosis', options.WeightKurtosis, ...
                'WeightVariance', options.WeightVariance, ...
                'WeightPulseNoise', options.WeightPulseNoise, ...
                'WeightPNR', options.WeightPNR, ...
                'WeightRatePenalty', options.WeightRatePenalty, ...
                'MaxNonPenalizedRate', options.MaxNonPenalizedRate);

        update_ipt_axes_graphics();
        fig.Pointer = "arrow";
    end

    function update_ipt_axes_graphics(~,~)
        fig.UserData.IPTLine.YData = fig.UserData.IPTs{fig.UserData.iRep}(fig.UserData.iPT(fig.UserData.iRep),:);
        xx = fig.UserData.t(fig.UserData.MUPulses{fig.UserData.iRep}{fig.UserData.iPT(fig.UserData.iRep)});
        yy = fig.UserData.IPTs{fig.UserData.iRep}(fig.UserData.iPT(fig.UserData.iRep),fig.UserData.MUPulses{fig.UserData.iRep}{fig.UserData.iPT(fig.UserData.iRep)});
        set(fig.UserData.MUPulseLine, 'XData', xx, 'YData', yy);
        set(fig.UserData.IPTThresholdLine,'Value',fig.UserData.IPTThreshold{fig.UserData.iRep}(fig.UserData.iPT(fig.UserData.iRep)));
        set(fig.UserData.Title, 'String', sprintf("IPT-%d | Rep-%d | N = %d | PNR = %4.1f dB", ...
            fig.UserData.iPT(fig.UserData.iRep), fig.UserData.iRep, numel(yy), round(fig.UserData.PNR{fig.UserData.iRep}(fig.UserData.iPT(fig.UserData.iRep)),1)));
    end

    function ipt_threshold_line_grab_callback(~, ~)
        fig.WindowButtonUpFcn = @ipt_threshold_line_release_callback;
        fig.UserData.IPTThresholdLine.ButtonDownFcn = [];
    end

    function ipt_threshold_line_release_callback(~, ~)
        fig.WindowButtonUpFcn = [];
        newThreshold = fig.UserData.IPTAxes.CurrentPoint(1,2);
        fig.UserData.IPTThresholdLine.Value = newThreshold;
        fig.UserData.IPTThreshold{fig.UserData.iRep}(fig.UserData.iPT(fig.UserData.iRep)) = newThreshold;
        fig.Pointer = "watch";
        drawnow();
        recalculate_ipts();
        fig.Pointer = "arrow";
        fig.UserData.IPTThresholdLine.ButtonDownFcn = @ipt_threshold_line_grab_callback;
    end

    function handle_changing_ipt_index(src,~)
        fig.UserData.iPT(fig.UserData.iRep) = src.Value;
        update_ipt_axes_graphics();
    end

    function handle_adding_new_rep(~,~)
        u_hat = fig.UserData.uni{end};
        p_hat = fig.UserData.proj{end};
        c_hat = fig.UserData.coeff{end};
        fig.UserData.iPT(end+1) = 1;
        fig.UserData.nRep = fig.UserData.nRep + 1;
        fig.UserData.iRep = fig.UserData.iRep + 1;
        fig.UserData.IPTThreshold{end+1} = fig.UserData.IPTThreshold{end};
        fig.UserData.Cost{end+1} = fig.UserData.Cost{end};
        fig.UserData.PNR{end+1} = fig.UserData.PNR{end};
        fig.UserData.CoV{end+1} = fig.UserData.CoV{end};
        fig.UserData.L2{end+1} = fig.UserData.L2{end};
        fig.UserData.Kurtosis{end+1} = fig.UserData.Kurtosis{end};
        fig.UserData.Variance{end+1} = fig.UserData.Variance{end};
        fig.UserData.PulseNoise{end+1} = fig.UserData.PulseNoise{end};
        fig.UserData.RatePenalty{end+1} = fig.UserData.RatePenalty{end};
        [fig.UserData.coeff{end+1}, fig.UserData.proj{end+1}] = ckc.uni_2_pcs(u_hat);
        fig.UserData.nPT(end+1) = size(fig.UserData.coeff{end},2);
        
        
        [fig.UserData.IPTs{end}, fig.UserData.MUPulses{end}, fig.UserData.uni{end+1}] = ...
                ckc.compute_IPT_fast(u_hat, p_hat, c_hat, ...
                    'SigmaThreshold', options.SigmaThreshold, ...
                    'IPTThreshold', fig.UserData.IPTThreshold{end}, ...
                    'MaxThresholdAdjustments', options.MaxThresholdAdjustments, ...
                    'ThresholdAdjustmentScalar', options.ThresholdAdjustmentScalar);
        fig.UserData.IPTs{end+1} = fig.UserData.IPTs{end};
        fig.UserData.MUPulses{end+1} = fig.UserData.MUPulses{end};

        recalculate_ipts();
    end

end