function PNR = compute_PNR(MUPulses,IPT,fsamp)
%COMPUTE_PNR  Calculates the pulse-to-noise ratio (PNR) 

if ~isempty(MUPulses)
    
    tmpMUPulses = [];
    for t = -round(3/2048*fsamp):round(3/2048*fsamp)
        tmpMUPulses = [tmpMUPulses MUPulses+t];
    end
    
    nInd = setdiff(MUPulses(1):MUPulses(end),tmpMUPulses);
    IPT = IPT/mean(IPT(MUPulses));
    
    tmpn = IPT(nInd);
    tmpn = tmpn(~isnan(tmpn));
    tmpn = tmpn(tmpn>=0);
    PNR = round(10*10*log10(  mean(IPT(MUPulses).^2) / mean(tmpn.^2)))/10;
else
    PNR = 0;
end
