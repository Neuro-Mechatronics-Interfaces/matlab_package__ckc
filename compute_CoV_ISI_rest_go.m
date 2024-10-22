function [CoV_ISI_Rest, ISI_Rest, CoV_ISI_Go, ISI_Go] = compute_CoV_ISI_rest_go(MUPulses, sync, fs, minSamplesInEpoch)

if nargin < 3
    fs = 2000;
end

if nargin < 4
    minSamplesInEpoch = 4000;
end

HIGH = find(sync);
LOW = find(~sync);
deltaHIGH = diff(HIGH);
deltaLOW = diff(LOW);
iGoRising = HIGH([false,deltaHIGH>1]);
iGoFalling = LOW([false,deltaLOW>1]);

if iGoFalling(1) < iGoRising(1)
    iGoFalling(1) = [];
end

if iGoFalling(end) < iGoRising(end)
    iGoFalling(end+1) = numel(sync);
    iGoFalling = circshift(iGoFalling,-1);
end

iRestFalling = iGoRising(2:end)-1;
iRestRising = iGoFalling(2:end)+1;

if ~iscell(MUPulses)
    MUPulses = {MUPulses};
end
CoV_ISI_Rest = nan(numel(MUPulses),1);
CoV_ISI_Go = nan(numel(MUPulses),1);
ISI_Rest = cell(numel(MUPulses),1);
ISI_Go = cell(numel(MUPulses),1);
for ii = 1:numel(MUPulses)
    ISI_Go{ii} = [];
    ISI_Rest{ii} = [];
    for ik = 1:numel(iGoRising)
        if iGoFalling(ik)-iGoRising(ik)<minSamplesInEpoch
            continue;
        end
        goPulses = MUPulses{ii}(ismember(MUPulses{ii},iGoRising(ik):iGoFalling(ik)));
        ISI_Go{ii} = [ISI_Go{ii}, diff(goPulses./fs)];
    end
    CoV_ISI_Go(ii) = std(ISI_Go{ii})/mean(ISI_Go{ii});
    for ik = 1:numel(iRestRising)
        if iRestFalling(ik)-iRestRising(ik)<minSamplesInEpoch
            continue;
        end
        restPulses = MUPulses{ii}(ismember(MUPulses{ii},iRestRising(ik):iRestFalling(ik)));
        ISI_Rest{ii} = [ISI_Rest{ii}, diff(restPulses./fs)];
    end
    CoV_ISI_Rest(ii) = std(ISI_Rest{ii})/mean(ISI_Rest{ii});
end


end