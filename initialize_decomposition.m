function [IPTs, MUPulses, info, t, R_inv] = initialize_decomposition(uni, fsamp, options)
%INITIALIZE_DECOMPOSITION  Initialize the decomposition for CKC
%
% Example:
%   data = TMSiSAGA.Poly5.read('Max_2024_03_30_B_22.poly5');
%   uni = ckc.preproc__hpf_exclude_interp_del2(data.samples(2:65,:));   
%   [IPTs, MUPulses, info, t, R_inv] = ckc.initialize_decomposition(uni, data.sample_rate);
%
% Syntax:
%   [IPTs, MUPulses, info, t, R_inv] = ckc.initialize_decomposition(uni, fsamp, 'Name', value, ...);
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
    options.ExtensionFactor = 40;
    options.NumPCs (1,:) {mustBePositive, mustBeInteger} = 7:3:10;
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
    options.MaxNonPenalizedRate (1,1) double = 60;
    options.MaxStrikes (1,1) {mustBeInteger, mustBePositive} = 10;
    options.PreferredIntervalTolerance (1,1) double = 0.2;
    options.InverseCovarianceMatrix = [];
    options.Verbose (1,1) logical = true;
end

warning('off','signal:findpeaks:largeMinPeakHeight');
if ~any(isnan(options.SubSampleVector))
    Y = uni(:,options.SubSampleVector(1):options.SubSampleVector(2));
    if isempty(options.InverseCovarianceMatrix)
        Ye = ckc.extend(Y,options.ExtensionFactor);
        Ry = Ye*Ye';
        R_inv = pinv(Ry);
    else
        R_inv = options.InverseCovarianceMatrix;
    end
    Ue = ckc.extend(uni,options.ExtensionFactor);
else
    Ue = ckc.extend(uni,options.ExtensionFactor);
    if isempty(options.InverseCovarianceMatrix)
        Ru = Ue*Ue';
        R_inv = pinv(Ru);
    else
        R_inv = options.InverseCovarianceMatrix;
    end
end

t = (0:(size(Ue,2)-1))/fsamp;

info = cell(1,numel(options.NumPCs));
MUPulses = cell(1,numel(options.NumPCs));
IPTs = cell(1,numel(options.NumPCs));
cost = cell(1,numel(options.NumPCs));
for iSet = 1:numel(options.NumPCs)
    [~, all_locs, explained] = ckc.uni_2_pks(uni, 'NPCs', options.NumPCs(iSet), 'Verbose', options.Verbose);
    IPTs{iSet} = nan(options.NumPCs(iSet), size(Ue,2));
    MUPulses{iSet} = cell(1, options.NumPCs(iSet));
    IPTThreshold = ones(1,options.NumPCs(iSet)).*options.IPTThreshold;
    for ii = 1:options.NumPCs(iSet)
        [IPTs{iSet}(ii,:), MUPulses{iSet}{ii}, IPTThreshold(ii)] = ckc.compute_IPT(Ue, R_inv, all_locs{ii}, ...
            'SigmaThreshold', options.SigmaThreshold, ...
            'ExtensionFactor', options.ExtensionFactor, ...
            'IPTThreshold', IPTThreshold(ii), ...
            'MaxThresholdAdjustments', options.MaxThresholdAdjustments, ...
            'ThresholdAdjustmentScalar', options.ThresholdAdjustmentScalar);    
    end
    [cost{iSet}, PNR, CoV, L2, Kurtosis, Variance, PulseNoise, RatePenalty] = ckc.compute_decomposition_cost(...
        uni, MUPulses{iSet}, IPTs{iSet}, fsamp, ...
        'WeightL2Norm', options.WeightL2Norm, ...
        'NoiseUpperBound', options.IPTNoiseUpperBound, ...
        'WeightKurtosis', options.WeightKurtosis, ...
        'WeightVariance', options.WeightVariance, ...
        'WeightPulseNoise', options.WeightPulseNoise, ...
        'WeightCoV', options.WeightCoV, ...
        'WeightPNR', options.WeightPNR, ...
        'WeightRatePenalty', options.WeightRatePenalty, ...
        'MaxNonPenalizedRate', options.MaxNonPenalizedRate);
    
    IPTConverged = false(1,options.NumPCs(iSet));
    IPTWorse = false(1,options.NumPCs(iSet));
    NStrikes = zeros(1,options.NumPCs(iSet));
    
    if options.Verbose
        fprintf(1,'Running IPT cleaning iterations...000%%\n');
    end
    
    for iIter = 1:options.IPTCleaningIterationsMax
        for ii = 1:options.NumPCs(iSet)
            if IPTConverged(ii)
                continue;
            end
            [IPTnew, newPulses,IPTThreshold(ii)] = ckc.compute_IPT(Ue, R_inv, MUPulses{iSet}{ii}, ...
                    'IPTThreshold', IPTThreshold(ii), ...
                    'ExtensionFactor', options.ExtensionFactor, ...
                    'SigmaThreshold', options.SigmaThreshold, ...
                    'MaxThresholdAdjustments', options.MaxThresholdAdjustments, ...
                    'ThresholdAdjustmentScalar', options.ThresholdAdjustmentScalar);
            [newCost, newPNR, newCoV, newL2, newKurtosis, newVariance, newPulseNoise, newRatePenalty] = ckc.compute_decomposition_cost(...
                uni, {newPulses}, IPTnew, fsamp, ...
                    'NoiseUpperBound', options.IPTNoiseUpperBound, ...
                    'WeightL2Norm', options.WeightL2Norm, ...
                    'WeightCoV', options.WeightCoV, ...
                    'WeightKurtosis', options.WeightKurtosis, ...
                    'WeightVariance', options.WeightVariance, ...
                    'WeightPulseNoise', options.WeightPulseNoise, ...
                    'WeightPNR', options.WeightPNR, ...
                    'WeightRatePenalty', options.WeightRatePenalty, ...
                    'MaxNonPenalizedRate', options.MaxNonPenalizedRate);
            IPTWorse(ii) = newCost > cost{iSet}(ii);
            if IPTWorse(ii)
                IPTThreshold(ii) = IPTThreshold(ii) - 1.5*options.IPTThresholdStepSize;
                
                Delta = ckc.compute_nearest_deltas(MUPulses{iSet}{ii},2);
                MedDelta = median(Delta,2);
                dMedDelta = MedDelta(2) - MedDelta(1);
                iLow = find(Delta(1,:) < (dMedDelta - options.PreferredIntervalTolerance*dMedDelta));
                if ~isempty(iLow)
                    iRemove = false(size(iLow));
                    for ik = 1:numel(iLow)
                        if abs(Delta(2,iLow(ik)) - MedDelta) < (dMedDelta * options.PreferredIntervalTolerance)
                            iRemove(ik) = true;
                        end
                    end  
                    MUPulses{iSet}{ii}(iLow(iRemove)) = [];
                else
                    [~,iMax] = max(IPTs{iSet}(ii,MUPulses{iSet}{ii}));
                    MUPulses{iSet}{ii}(iMax) = [];
                end
                
                NStrikes(ii) = NStrikes(ii) + 1;
                if NStrikes(ii) == options.MaxStrikes
                    IPTConverged(ii) = true;
                    if options.Verbose
                        fprintf(1,'\b\b\b\b\b\nIPT-%d converged (iter=%d).\n', ii, iIter);
                        fprintf(1,'Running...%03d%%\n', round(iIter*100 / options.IPTCleaningIterationsMax));
                    end
                end
            else
                cost{iSet}(ii) = newCost;
                PNR(ii) = newPNR;
                CoV(ii) = newCoV;
                L2(ii) = newL2;
                MUPulses{iSet}{ii} = newPulses;
                Kurtosis(ii) = newKurtosis;
                Variance(ii) = newVariance;
                PulseNoise(ii) = newPulseNoise;
                RatePenalty(ii) = newRatePenalty;
                IPTs{iSet}(ii,:) = IPTnew;
                IPTThreshold(ii) = IPTThreshold(ii) + options.IPTThresholdStepSize;
            end
        end
        if options.Verbose
            fprintf(1,'\b\b\b\b\b%03d%%\n', round(iIter*100 / options.IPTCleaningIterationsMax));
        end
    end
    info{iSet} = struct('num_pcs', options.NumPCs(iSet), 'explained', explained, 'cost', cost{iSet},  ...
        'CostTerms', struct2table(struct('PNR', PNR', 'CoV', CoV', 'L2', L2', 'Kurtosis', Kurtosis', 'Variance', Variance', 'PulseNoise', PulseNoise', 'RatePenalty', RatePenalty')), ...
        'WeightedCostTerms', struct2table(struct('PNR', PNR' * options.WeightPNR, 'CoV', options.WeightCoV * CoV', 'L2', options.WeightL2Norm * L2', 'Kurtosis', options.WeightKurtosis * Kurtosis', 'Variance', options.WeightVariance * Variance', 'PulseNoise', options.WeightPulseNoise * PulseNoise', 'RatePenalty', RatePenalty' * options.WeightRatePenalty)), ...
        'parameters', options);
end
warning('on','signal:findpeaks:largeMinPeakHeight');
bestVal = inf;
iBest = 0;
for ii = 1:numel(info)
    meanCost = mean(info{ii}.cost(~isnan(info{ii}.cost)));
    if meanCost < bestVal
        bestVal = meanCost;
        iBest = ii;
    end
end
info = info{iBest};
IPTs = IPTs{iBest};
MUPulses = MUPulses{iBest};
end