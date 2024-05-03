function [MUAP, MUAPrms, MUAPp2p, MUAPduration,AllMUAPs] = extractMUAPs(SIG,MUPulses,fsamp)
%EXTRACTMUAPS Function to actually get individual MUAP waveforms. 

AllMUAPs = cell(size(SIG));
MUAP = cell(size(SIG));
MUAPrms = nan(size(SIG));
MUAPp2p = nan(size(SIG));
MUAPduration = nan(size(SIG));
for r=1:size(SIG,1)
    for c=1:size(SIG,2)
        if ~isempty(SIG{r,c})
            
            tmpY = SIG{r,c}; % SD (DD when SIG is SD)
            tmpMUAP=ckc.extractMUAPocurrences(tmpY,MUPulses,round(0.03*fsamp));
            AllMUAPs{r,c} = tmpMUAP{1};
            MUAP{r,c} = mean(AllMUAPs{r,c});
            MUAPrms(r,c) = sqrt(mean(MUAP{r,c}.^2));
            MUAPp2p(r,c) = max(MUAP{r,c}) - min(MUAP{r,c});
            
            ind = find(abs(MUAP{r,c}) >= 0.8 * max(abs(MUAP{r,c})));
            MUAPduration(r,c) = max(ind) - min(ind);
        end
    end
end
